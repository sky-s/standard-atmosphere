function varargout = nsatmos(h, varargin)
%  NSATMOS  Find gas properties in non-standard atmospheres.
%   [rho,a,T,P,nu,z,sigma] = NSATMOS(h,varargin)
%
%   NSATMOS by itself gives atmospheric properties at sea level on a standard day.
%
%   NSATMOS(h) returns the properties of the 1976 Standard Atmosphere at
%   geopotential altitude h, where h is a scalar, vector, matrix, or ND array.
%   This is the same as ATMOS(h) for backward compatibility.
% 
%   The input h can be followed by parameter/value pairs for further control of
%   NSATMOS. Possible parameters are:
%     atmosphere   - Type of non-standard atmosphere. Options:
%                    'standard' (default) - 1976 Standard Atmosphere
%                    'hot'      - Hot day atmosphere (ISA +20°C to surface)
%                    'cold'     - Cold day atmosphere (ISA -20°C to surface)
%                    'tropical' - Tropical atmosphere profile
%                    'polar'    - Polar atmosphere profile
%     tOffset      - Additional temperature offset (same as ATMOS)
%     tAbsolute    - Absolute air temperature override (same as ATMOS)
%     altType      - Specify type of input altitude, either 'geopotential' (h)
%                    or 'geometric' (z). Default altType = 'geopotential'.
%     structOutput - When set, NSATMOS produces a single struct output with fields
%                    rho, a, T, P, nu, and either z or h (whichever complements
%                    input altType). Default structOutput = false.
%     units        - String for units of inputs and outputs, either 'SI'
%                    (default) or 'US'. This is ignored if the provided input h
%                    is a DimVar, in which case all outputs are also DimVars.
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
%   Example 1: Compare atmospheric properties at 10 km for different atmospheres
%       h = 10000;
%       [rho_std,a_std,T_std] = nsatmos(h, 'atmosphere', 'standard');
%       [rho_hot,a_hot,T_hot] = nsatmos(h, 'atmosphere', 'hot');
%       [rho_cold,a_cold,T_cold] = nsatmos(h, 'atmosphere', 'cold');
%       fprintf('At 10 km: Standard T=%.1f K, Hot T=%.1f K, Cold T=%.1f K\n', ...
%               T_std, T_hot, T_cold);
%
%   Example 2: Plot temperature profiles for different atmosphere types
%       h = 0:1000:20000;
%       T_std = nsatmos(h, 'atmosphere', 'standard', 'structOutput', true);
%       T_hot = nsatmos(h, 'atmosphere', 'hot', 'structOutput', true);
%       T_tropical = nsatmos(h, 'atmosphere', 'tropical', 'structOutput', true);
%       plot(T_std.T, h/1000, T_hot.T, h/1000, T_tropical.T, h/1000);
%       legend('Standard', 'Hot', 'Tropical');
%       xlabel('Temperature (K)'); ylabel('Altitude (km)');
%
%   References: 
%     MIL-STD-210C, MIL-HDBK-310, www.pdas.com/milstd210.html
%     U.S. Standard Atmosphere, 1976
%
%   See also ATMOS, ATMOSISA, TROPOS.
%
%   [rho,a,T,P,nu,z,sigma] = NSATMOS(h,varargin)

%   Copyright 2024 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715
% 
%   References: MIL-STD-210C; MIL-HDBK-310; www.pdas.com/milstd210.html

%% User-customizable defaults:
defaultUnits = 'SI'; % Alternate: 'US'
defaultStructOutput = false;
defaultAtmosphere = 'standard';

%% Parse inputs:
if nargin == 0
    h = 0;
end

% Quick return of sea level conditions for default case
if nargin <= 1 && ~nnz(h) && ischar(h) == false
    % For standard atmosphere at sea level
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

% Parse inputs
p = inputParser;
addParameter(p,'atmosphere',defaultAtmosphere);
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
else
    error('Invalid units. Expected: ''SI'' or ''US''.')
end    

if strcmpi(p.Results.altType,'geopotential')
    geomFlag = false;
elseif strcmpi(p.Results.altType,'geometric')
    geomFlag = true;
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
end
if isa(tAbsolute,'DimVar')
    tAbsolute = tAbsolute/u.K;
end

if convertUnits
    h = h * 0.3048;
    tOffset   = tOffset   * 5/9;
    tAbsolute = tAbsolute * 5/9;
end

%% Define atmosphere profiles:
% Standard atmosphere (same as 1976 Standard Atmosphere)
switch lower(atmosphere)
    case 'standard'
        %  Lapse rate Base Temp       Base Geop. Alt    Base Pressure
        %   Ki (°C/m) Ti (°K)         Hi (m)            P (Pa)
        D = [-0.0065   288.15          0                 101325            % Troposphere
             0         216.65          11000             22632.04059693474 % Tropopause
             0.001     216.65          20000             5474.877660660026 % Stratosph. 1
             0.0028    228.65          32000             868.0158377493657 % Stratosph. 2
             0         270.65          47000             110.9057845539146 % Stratopause
             -0.0028   270.65          51000             66.938535373039073% Mesosphere 1
             -0.002    214.65          71000             3.956392754582863 % Mesosphere 2
             0         186.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'hot'
        % Hot day atmosphere - Standard atmosphere with +20°C at surface
        % Temperature tapers from +20°C at surface to standard at tropopause
        D = [-0.0045   308.15          0                 101325            % Troposphere (reduced lapse)
             0         236.65          11000             22632.04059693474 % Tropopause (+20K)
             0.001     236.65          20000             5474.877660660026 % Stratosph. 1
             0.0028    248.65          32000             868.0158377493657 % Stratosph. 2
             0         290.65          47000             110.9057845539146 % Stratopause
             -0.0028   290.65          51000             66.938535373039073% Mesosphere 1
             -0.002    234.65          71000             3.956392754582863 % Mesosphere 2
             0         206.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'cold'
        % Cold day atmosphere - Standard atmosphere with -20°C at surface
        D = [-0.0085   268.15          0                 101325            % Troposphere (increased lapse)
             0         196.65          11000             22632.04059693474 % Tropopause (-20K)
             0.001     196.65          20000             5474.877660660026 % Stratosph. 1
             0.0028    208.65          32000             868.0158377493657 % Stratosph. 2
             0         250.65          47000             110.9057845539146 % Stratopause
             -0.0028   250.65          51000             66.938535373039073% Mesosphere 1
             -0.002    194.65          71000             3.956392754582863 % Mesosphere 2
             0         166.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'tropical'
        % Tropical atmosphere - Higher surface temperature, modified lapse
        % Based on MIL-STD-210C tropical profiles
        D = [-0.0055   303.15          0                 101325            % Troposphere (higher surface temp)
             0         241.65          11000             22632.04059693474 % Tropopause
             0.001     241.65          20000             5474.877660660026 % Stratosph. 1
             0.0028    253.65          32000             868.0158377493657 % Stratosph. 2
             0         295.65          47000             110.9057845539146 % Stratopause
             -0.0028   295.65          51000             66.938535373039073% Mesosphere 1
             -0.002    239.65          71000             3.956392754582863 % Mesosphere 2
             0         211.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
    case 'polar'
        % Polar atmosphere - Lower surface temperature, modified lapse
        % Based on MIL-STD-210C polar profiles
        D = [-0.0095   263.15          0                 101325            % Troposphere (lower surface temp)
             0         191.65          11000             22632.04059693474 % Tropopause
             0.001     191.65          20000             5474.877660660026 % Stratosph. 1
             0.0028    203.65          32000             868.0158377493657 % Stratosph. 2
             0         245.65          47000             110.9057845539146 % Stratopause
             -0.0028   245.65          51000             66.938535373039073% Mesosphere 1
             -0.002    189.65          71000             3.956392754582863 % Mesosphere 2
             0         161.94590831019 84852.04584490575 0.373377242877530];% Mesopause
        
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
T = D(:,2);	%°K
H = D(:,3);	%m
P = D(:,4);	%Pa

R = P(1)/T(1)/rho0; %N-m/kg-K
% Note: Using sea level density ratio, pressure may need adjustment for non-standard atmospheres
% For accuracy, recalculate based on actual surface conditions
R = 287.0528742470439; % Use standard gas constant

%% Convert from geometric altitude to geopotential altitude, if necessary.
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

%% Switch between using calculated temp and provided absolute temp.
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