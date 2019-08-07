function varargout = atmos(h,varargin)
%  ATMOS  Find gas properties in the 1976 Standard Atmosphere.
%   [rho,a,T,P,nu,z,sigma] = ATMOS(h,varargin)
%
%   ATMOS by itself gives atmospheric properties at sea level on a standard day.
%
%   ATMOS(h) returns the properties of the 1976 Standard Atmosphere at
%   geopotential altitude h, where h is a scalar, vector, matrix, or ND array.
% 
%   The input h can be followed by parameter/value pairs for further control of
%   ATMOS. Possible parameters are:
%     tOffset      - Returns properties when the temperature is tOffset degrees 
%                    above or below standand conditions. h and tOffset must be
%                    the same size or else one must be a scalar. Default is no
%                    offset. Note that this is an offset, so when converting
%                    between Celsius and Fahrenheit, use only the scaling factor
%                    (dC/dF = dK/dR = 5/9).
%     tAbsolute    - Similar to tOffest, but an absolute air temperature is 
%                    provided (°K or °R) instead of an offset from the standard 
%                    temperature. Supersedes tOffset if both are provided.
%     altType      - Specify type of input altitude, either 'geopotential' (h)
%                    or 'geometric' (z). Default altType = 'geopotential'.
%     structOutput - When set, ATMOS produces a single struct output with fields
%                    rho, a, T, P, nu, and either z or h (whichever complements
%                    input altType). Default structOutput = false.
%     units        - String for units of inputs and outpus, either 'SI'
%                    (default) or 'US'. This is ignored if the provided input h
%                    is a DimVar, in which case all outputs are also DimVars and
%                    expected tOffset is either a DimVar or in °C/°K.
%                                 Description:         SI:           US:
%                     Input:      --------------       -----         -----
%                       h | z     Altitude or height   m             ft
%                       tOffset   Temp. offset         °C/°K         °F/°R
%                     Output:     --------------       -----         -----
%                       rho       Density              kg/m^3        slug/ft^3
%                       a         Speed of sound       m/s           ft/s
%                       T         Temperature          °K            °R
%                       P         Pressure             Pa            lbf/ft^2
%                       nu        Kinem. viscosity     m^2/s         ft^2/s
%                       z | h     Height or altitude   m             ft
%                       sigma     Density ratio        -             -
%
%   ATMOS returns properties the same size as h and/or tOffset (P does not vary
%   with temperature offset and is always the size of h).
%
%   Example 1: Find atmospheric properties at every 100 m of geometric height
%   for an off-standard atmosphere with temperature offset varying +/- 25°C
%   sinusoidally with a period of 4 km.
%       z = 0:100:86000;
%       [rho,a,T,P,nu,h,sigma] = atmos(z,'tOffset',25*sin(pi*z/2000),...
%                                        'altType','geometric');
%       semilogx(sigma,h/1000)
%       title('Density variation with sinusoidal off-standard atmosphere')
%       xlabel('Density ratio, \sigma'); ylabel('Geopotential altitude (km)')
%
%   Example 2: Create tables of atmospheric properties up to 30,000 ft for a
%   cold (-20°C), standard, and hot (+20°C) day with columns
%   [h(ft) z(ft) rho(slug/ft³) sigma a(ft/s) T(R) P(psf) µ(slug/ft-s) nu(ft²/s)]
%   leveraging n-dimensional array capability.
%       [~,h,dT] = meshgrid(0,-5000:1000:30000,[-20 0 20]);
%       [rho,a,T,P,nu,z,sigma] = atmos(h,'tOffset',dT*9/5,'units','US');
%       t = [h z rho sigma a T P nu.*rho nu];
%       format short e
%       varNames = {'h' 'z' 'rho' 'sigma' 'a' 'T' 'P' 'mu' 'nu'};
%       ColdTable       = array2table(t(:,:,1),'VariableNames',varNames)
%       StandardTable   = array2table(t(:,:,2),'VariableNames',varNames)
%       HotTable        = array2table(t(:,:,3),'VariableNames',varNames)
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
%   above 86 km (with a lapse rate of 0°/km) and below 0.
%
%   See also ATMOSISA, ATMOSNONSTD, TROPOS,
%     DENSITYALT - http://www.mathworks.com/matlabcentral/fileexchange/39325,
%     UNITS      - http://www.mathworks.com/matlabcentral/fileexchange/38977.
%
%   [rho,a,T,P,nu,z,sigma] = ATMOS(h,varargin)

%   Copyright 2015 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715
% 
%   References: ESDU 77022; www.pdas.com/atmos.html

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
    a = 340.2940;
    temp = 288.15;
    press = 101325;
    kvisc = 1.4607186e-05;
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
addParameter(p,'tOffset',0);
addParameter(p,'tAbsolute',[]);
addParameter(p,'units',defaultUnits);
addParameter(p,'altType','geopotential');
addParameter(p,'structOutput',defaultStructOutput);
parse(p,varargin{:});

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


%% Constants, etc.:

%  Lapse rate Base Temp       Base Geop. Alt    Base Pressure
%   Ki (°C/m) Ti (°K)         Hi (m)            P (Pa)
D =[-0.0065   288.15          0                 101325            % Troposphere
    0         216.65          11000             22632.04059693474 % Tropopause
    0.001     216.65          20000             5474.877660660026 % Stratosph. 1
    0.0028    228.65          32000             868.0158377493657 % Stratosph. 2
    0         270.65          47000             110.9057845539146 % Stratopause
    -0.0028   270.65          51000             66.938535373039073% Mesosphere 1
    -0.002    214.65          71000             3.956392754582863 % Mesosphere 2
    0         186.94590831019 84852.04584490575 .373377242877530];% Mesopause

% Constants:
rho0 = 1.225;   % Sea level density, kg/m^3
gamma = 1.4;
g0 = 9.80665;   %m/sec^2
RE = 6356766;   %Radius of the Earth, m
Bs = 1.458e-6;  %N-s/m2 K1/2
S = 110.4;      %K

K = D(:,1);	%°K/m
T = D(:,2);	%°K
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
    
    
    if K(i) == 0 % No temperature lapse.
        temp(n) = T(i);
        press(n) = P(i) * exp(-g0*(hGeop(n)-H(i))/(T(i)*R));
    else
        TonTi = 1 + K(i)*(hGeop(n) - H(i))/T(i);
        temp(n) = TonTi*T(i); 
        press(n) = P(i) * TonTi.^(-g0/(K(i)*R)); % Undefined for K = 0.
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