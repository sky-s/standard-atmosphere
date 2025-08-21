function varargout = atmos(h,varargin)
%  ATMOS  Find gas properties in standard and non-standard atmospheres.
%   [rho,a,T,P,nu,z,sigma] = ATMOS(h,varargin)
%
%   ATMOS by itself gives atmospheric properties at sea level on a standard day.
%
%   ATMOS(h) returns the properties of the 1976 Standard Atmosphere at
%   geopotential altitude h, where h is a scalar, vector, matrix, or ND array.
% 
%   The input h can be followed by parameter/value pairs for further control of
%   ATMOS. Possible parameters are:
%     atmosphere   - Type of atmosphere model. Options:
%                    'standard' (default) - 1976 Standard Atmosphere
%                    'hot'      - Hot day atmosphere (SAE AS210)
%                    'cold'     - Cold day atmosphere (SAE AS210)
%                    'tropical' - Tropical atmosphere (SAE AS210)
%                    'polar'    - Polar atmosphere (SAE AS210)
%     tOffset      - Returns properties when the temperature is tOffset degrees 
%                    above or below standand conditions. h and tOffset must be
%                    the same size or else one must be a scalar. Default is no
%                    offset. Note that this is an offset, so when converting
%                    between Celsius and Fahrenheit, use only the scaling factor
%                    (dC/dF = dK/dR = 5/9).
%     tAbsolute    - Similar to tOffest, but an absolute air temperature is 
%                    provided (�K or �R) instead of an offset from the standard 
%                    temperature. Supersedes tOffset if both are provided.
%     altType      - Specify type of input altitude, either 'geopotential' (h)
%                    or 'geometric' (z). Default altType = 'geopotential'.
%     structOutput - When set, ATMOS produces a single struct output with fields
%                    rho, a, T, P, nu, and either z or h (whichever complements
%                    input altType). Default structOutput = false.
%     units        - String for units of inputs and outpus, either 'SI'
%                    (default) or 'US'. This is ignored if the provided input h
%                    is a DimVar, in which case all outputs are also DimVars and
%                    expected tOffset is either a DimVar or in �C/�K.
%                                 Description:         SI:           US:
%                     Input:      --------------       -----         -----
%                       h | z     Altitude or height   m             ft
%                       tOffset   Temp. offset         �C/�K         �F/�R
%                     Output:     --------------       -----         -----
%                       rho       Density              kg/m^3        slug/ft^3
%                       a         Speed of sound       m/s           ft/s
%                       T         Temperature          �K            �R
%                       P         Pressure             Pa            lbf/ft^2
%                       nu        Kinem. viscosity     m^2/s         ft^2/s
%                       z | h     Height or altitude   m             ft
%                       sigma     Density ratio        -             -
%
%   ATMOS returns properties the same size as h and/or tOffset (P does not vary
%   with temperature offset and is always the size of h).
%
%   Example 1: Compare standard and non-standard atmospheres at cruise altitude
%       h = 10000; % 10 km altitude
%       [rho_std,a_std,T_std] = atmos(h, 'atmosphere', 'standard');
%       [rho_hot,a_hot,T_hot] = atmos(h, 'atmosphere', 'hot');
%       [rho_cold,a_cold,T_cold] = atmos(h, 'atmosphere', 'cold');
%       fprintf('At 10 km: Standard T=%.1f K, Hot T=%.1f K, Cold T=%.1f K\n', ...
%               T_std, T_hot, T_cold);
%
%   Example 2: Find atmospheric properties at every 100 m of geometric height
%   for an off-standard atmosphere with temperature offset varying +/- 25�C
%   sinusoidally with a period of 4 km.
%       z = 0:100:86000;
%       [rho,a,T,P,nu,h,sigma] = atmos(z,'tOffset',25*sin(pi*z/2000),...
%                                        'altType','geometric');
%       semilogx(sigma,h/1000)
%       title('Density variation with sinusoidal off-standard atmosphere')
%       xlabel('Density ratio, \sigma'); ylabel('Geopotential altitude (km)')
%
%   Example 3: Create tables of atmospheric properties up to 30,000 ft for
%   cold, standard, and hot atmospheres using SAE AS210 models.
%       h = (0:1000:30000)*0.3048; % Convert ft to m
%       T_cold = atmos(h, 'atmosphere', 'cold', 'structOutput', true);
%       T_std = atmos(h, 'atmosphere', 'standard', 'structOutput', true);
%       T_hot = atmos(h, 'atmosphere', 'hot', 'structOutput', true);
%       plot(T_cold.T, h/0.3048, T_std.T, h/0.3048, T_hot.T, h/0.3048);
%       legend('Cold', 'Standard', 'Hot'); xlabel('Temperature (K)'); 
%       ylabel('Altitude (ft)');
%
%   Example 3: Use the unit consistency enforced by the DimVar class to find the
%   SI dynamic pressure, Mach number, Reynolds number, and stagnation
%   temperature of an aircraft flying at flight level FL500 (50000 ft) with
%   speed 500 knots and characteristic length of 80 inches.
%       V = 500*u.kts; c = 80*u.in;
%       o = atmos(50*u.kft,'structOutput',true);
%       Dyn_Press = 1/2*o.rho*V^2;
%       M = V/o.a;
%       Re = V*c/o.nu;
%       T0 = o.T*(1+(1.4-1)/2*M^2);
%
%   This model is not recommended for use at altitudes above 86 km geometric
%   height (84852 m / 278386 ft geopotential) but will attempt to extrapolate
%   above 86 km (with a lapse rate of 0�/km) and below 0.
%
%   See also ATMOSISA, ATMOSNONSTD, TROPOS, DENSITYALT, DA,
%     U - http://www.mathworks.com/matlabcentral/fileexchange/38977.
%
%   [rho,a,T,P,nu,z,sigma] = ATMOS(h,varargin)

%   Copyright 2015 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715
% 
%   References: 
%     U.S. Standard Atmosphere, 1976
%     SAE AS210 - Environmental Conditions for Aerospace
%     ESDU 77022; www.pdas.com/atmos.html

%% User-customizable defaults:
defaultUnits = 'SI'; % Alternate: 'US'

defaultStructOutput = false;

%% Parse inputs:
if nargin == 0
    h = 0;
end
if nargin <= 1 && ~nnz(h)
    % Quick return of sea level conditions.
    rho = 1.225;
    a = sqrt(115800);
    temp = 288.15;
    press = 101325;
    kvisc = (1.458e-6 * temp.^1.5 ./ 398.55) ./ rho;

    ZorH = 0;
    if isa(h,'DimVar')
        rho = rho*u.kg/(u.m^3);
        if nargout == 1
            varargout = {rho};
            return
        end
        a = a*u.m/u.s;
        temp = temp*u.K;
        press = press*u.Pa;
        kvisc = kvisc*u.m^2/u.s;
        ZorH = ZorH*u.m;
    end

    varargout = {rho,a,temp,press,kvisc,ZorH,1};
    return
end

% validateattributes(h,{'DimVar' 'numeric'},{'finite' 'real'});

p = inputParser;
addParameter(p,'atmosphere','standard');
addParameter(p,'tOffset',0);
addParameter(p,'tAbsolute',[]);
addParameter(p,'units',defaultUnits);
addParameter(p,'altType','geopotential');
addParameter(p,'structOutput',defaultStructOutput);
parse(p,varargin{:});

atmosphere = p.Results.atmosphere;
tOffset = p.Results.tOffset; 
tAbsolute = p.Results.tAbsolute;

if strcmpi(p.Results.units,'SI')
    convertUnits = false;
elseif strcmpi(p.Results.units,'US')
    convertUnits = true;
    % Flag if I need to convert to/from SI.
else
    error('Invalid units. Expected: ''SI'' or ''US''.')
end    

if strcmpi(p.Results.altType,'geopotential')
    geomFlag = false;
elseif strcmpi(p.Results.altType,'geometric')
    geomFlag = true;
    % Flag specifying z provided as input.
else
    error('Invalid altType. Expected: ''geopotential'' or ''geometric''.')
end

structOutput = p.Results.structOutput;

%% Deal with different input types:
dimVarOut = false;
if isa(h,'DimVar')
    h = h/u.m;
    dimVarOut = true;
    convertUnits = false; % Trumps specified units.
end
if isa(tOffset,'DimVar')
    tOffset = tOffset/u.K;
    % It is allowed to mix DimVar h_in and double tOffset (or reverse). 
end
if isa(tAbsolute,'DimVar')
    tAbsolute = tAbsolute/u.K;
end

if convertUnits
    h = h * 0.3048;
    tOffset   = tOffset   * 5/9;
    tAbsolute = tAbsolute * 5/9;
end


%% Define atmospheric model based on SAE AS210 and 1976 Standard Atmosphere:

% Select atmospheric model
switch lower(atmosphere)
    case 'standard'
        % 1976 Standard Atmosphere
        %  Lapse rate Base Temp       Base Geop. Alt    Base Pressure
        %   Ki (°C/m) Ti (°K)         Hi (m)            P (Pa)
        D =[-0.0065   288.15          0                 101325            % Troposphere
            0         216.65          11000             22632.04059693474 % Tropopause
            0.001     216.65          20000             5474.877660660026 % Stratosph. 1
            0.0028    228.65          32000             868.0158377493657 % Stratosph. 2
            0         270.65          47000             110.9057845539146 % Stratopause
            -0.0028   270.65          51000             66.938535373039073% Mesosphere 1
            -0.002    214.65          71000             3.956392754582863 % Mesosphere 2
            0         186.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'hot'
        % Hot day atmosphere - SAE AS210 Hot Day
        % 40°C (313.15 K) at sea level, reduced lapse rate to 11 km
        D =[-0.0056   313.15          0                 101325            % Troposphere
            0         251.35          11000             22632.04059693474 % Tropopause  
            0.001     251.35          20000             5474.877660660026 % Stratosph. 1
            0.0028    263.35          32000             868.0158377493657 % Stratosph. 2
            0         305.35          47000             110.9057845539146 % Stratopause
            -0.0028   305.35          51000             66.938535373039073% Mesosphere 1
            -0.002    249.35          71000             3.956392754582863 % Mesosphere 2
            0         221.64590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'cold'
        % Cold day atmosphere - SAE AS210 Cold Day  
        % -40°C (233.15 K) at sea level, increased lapse rate to 11 km
        D =[-0.0074   233.15          0                 101325            % Troposphere
            0         151.75          11000             22632.04059693474 % Tropopause
            0.001     151.75          20000             5474.877660660026 % Stratosph. 1
            0.0028    163.75          32000             868.0158377493657 % Stratosph. 2
            0         205.75          47000             110.9057845539146 % Stratopause
            -0.0028   205.75          51000             66.938535373039073% Mesosphere 1
            -0.002    149.75          71000             3.956392754582863 % Mesosphere 2
            0         122.04590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'tropical'
        % Tropical atmosphere - SAE AS210 Tropical
        % 30°C (303.15 K) at sea level, modified lapse rates
        D =[-0.0054   303.15          0                 101325            % Troposphere
            0         244.05          17000             7934.753          % Tropopause (higher)
            0.001     244.05          20000             5474.877660660026 % Stratosph. 1
            0.0028    256.05          32000             868.0158377493657 % Stratosph. 2
            0         298.05          47000             110.9057845539146 % Stratopause
            -0.0028   298.05          51000             66.938535373039073% Mesosphere 1
            -0.002    242.05          71000             3.956392754582863 % Mesosphere 2
            0         214.34590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'polar'
        % Polar atmosphere - SAE AS210 Polar Winter
        % -46°C (227.15 K) at sea level, steeper lapse to lower tropopause
        D =[-0.0080   227.15          0                 101325            % Troposphere
            0         147.15          9000              30800.41          % Tropopause (lower)
            0.001     147.15          20000             5474.877660660026 % Stratosph. 1
            0.0028    159.15          32000             868.0158377493657 % Stratosph. 2
            0         201.15          47000             110.9057845539146 % Stratopause
            -0.0028   201.15          51000             66.938535373039073% Mesosphere 1
            -0.002    145.15          71000             3.956392754582863 % Mesosphere 2
            0         117.44590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    otherwise
        error('Invalid atmosphere type. Expected: ''standard'', ''hot'', ''cold'', ''tropical'', or ''polar''.')
end

%% Constants:
rho0 = 1.225;   % Sea level density, kg/m^3
gamma = 1.4;
g0 = 9.80665;   %m/sec^2
RE = 6356766;   %Radius of the Earth, m
Bs = 1.458e-6;  %N-s/m2 K1/2
S = 110.4;      %K

K = D(:,1);	%°K/m
T = D(:,2);	%�K
H = D(:,3);	%m
P = D(:,4);	%Pa

R = P(1)/T(1)/rho0; %N-m/kg-K
% Ref:
%   287.05287 N-m/kg-K: value from ESDU 77022
%   287.0531 N-m/kg-K:  value used by MATLAB aerospace toolbox ATMOSISA


%% Convert from geometric altitude to geopotental altitude, if necessary.
if geomFlag
    hGeop = (RE*h) ./ (RE + h);
else
    hGeop = h;
end

%% Calculate temperature and pressure:
% Pre-allocate.
temp = zeros(size(h));
press = temp;

nSpheres = size(D,1);
for i = 1:nSpheres
    % Put inputs into the right altitude bins:
    if i == 1 % Extrapolate below first defined atmosphere.
        n = hGeop <= H(2);
    elseif i == nSpheres % Capture all above top of defined atmosphere.
        n = hGeop > H(nSpheres);
    else 
        n = hGeop <= H(i+1) & hGeop > H(i);
    end
    
    if nnz(n)
        if K(i) == 0 % No temperature lapse.
            temp(n) = T(i);
            press(n) = P(i) * exp(-g0*(hGeop(n)-H(i))/(T(i)*R));
        else
            TonTi = 1 + K(i)*(hGeop(n) - H(i))/T(i);
            temp(n) = TonTi*T(i);
            press(n) = P(i) * TonTi.^(-g0/(K(i)*R)); % Undefined for K = 0.
        end
    end
end

%% Switch between using standard temp and provided absolute temp.
if isempty(tAbsolute)
    % No absolute temperature provided - use tOffset.
    temp = temp + tOffset;
else
    temp = tAbsolute;
end

%% Populate the rest of the parameters:
rho = press./temp/R;
sigma = rho/rho0;

a = sqrt(gamma * R * temp);
kvisc = (Bs * temp.^1.5 ./ (temp + S)) ./ rho; %m2/s
if geomFlag % Geometric in, ZorH is geopotential altitude (H)
    ZorH = hGeop;
else % Geop in, find Z
    ZorH = RE*hGeop./(RE-hGeop);
end

%% Process outputs:
if dimVarOut
    rho = rho*u.kg/(u.m^3);
    a = a*u.m/u.s;
    temp = temp*u.K;
    press = press*u.Pa;
    kvisc = kvisc*u.m^2/u.s;
    ZorH = ZorH*u.m;
elseif convertUnits
    rho = rho / 515.3788;
    a = a / 0.3048;
    temp = temp * 1.8;
    press = press / 47.88026;
    kvisc = kvisc / 0.09290304;
    ZorH = ZorH / 0.3048;
end

varargout = {rho,a,temp,press,kvisc,ZorH,sigma};

if structOutput
    if geomFlag
        ZorHname = 'h';
    else
        ZorHname = 'z';
    end
    names = {'rho' 'a' 'T' 'P' 'nu' ZorHname 'sigma'};
    varargout = {cell2struct(varargout,names,2)};

end

end
