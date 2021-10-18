function h = da(rho,opts)
% DA  Returns geopotential altitude corresponding to the given array of air
% densities in the 1976 Standard Atmosphere.
%
%   h = DA(rho) returns geopotential altitude, h, in meters as a function of air
%   density, rho, provided in kg/m³.
% 
%   h = DA(rho,'units','US') returns h in feet for a given rho in slug/ft³. 
% 
%   If the input rho is a DimVar, specified units are ignored, and the output h
%   will also be a DimVar.
% 
%   DA is valid for the entire standard atmosphere up through the mesopause (86
%   km height). It assumes that all that is known is air density. If pressure or
%   temperature and density are known, there exist more straightforward methods
%   for calculating density altitude:
%       P = rho*R*T; h = h0 * (1 - P^0.190284), where P is in atmospheres and 
%           h0 = 145366.45 ft or 44307.694 m.
%       (http://www.srh.noaa.gov/images/epz/wxcalc/pressureAltitude.pdf)
% 
% 
%   See also ATMOS, DA, 
%     U - http://www.mathworks.com/matlabcentral/fileexchange/38977.

%   Copyright Sky Sartorius
%   Author contact: mathworks.com/matlabcentral/fileexchange/authors/101715

%   Acknowledgement: Rob McDonald for prompting this closed-form version of
%   densityalt.

arguments
    rho
    opts.units {mustBeMember(opts.units,{'SI','US'})} = 'SI'
end

if strcmpi(opts.units,'SI')
    convertUnits = false;
elseif strcmpi(opts.units,'US')
    convertUnits = true;
    % Flag if I need to convert to/from SI.
else
    error('Invalid units. Expected: ''SI'' or ''US''.')
end    

%% Deal with different input types:
dimVarOut = false;
if isa(rho,'DimVar')
    rho = rho*u.m3/u.kg;
    dimVarOut = true;
    convertUnits = false; % Trumps specified units.
end

if convertUnits
    rho = rho * 515.3788183931961;
end

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

% Constants:
rho0 = 1.225; % Sea level density, kg/m^3
g0 = 9.80665;   %m/sec^2

K = D(:, 1); %°K/m
T = D(:, 2); %°K
H = D(:, 3); %m
P = D(:, 4); %Pa

R = P(1) / T(1) / rho0; %N-m/kg-K

% Density at base of each atmosphere layer.
RHO = P ./ (T * R);

h = zeros(size(rho));

nSpheres = size(D,1);
for iSphere = 1:nSpheres
    % Put inputs into the right altitude bins:
    if iSphere == 1 % Extrapolate below first defined atmosphere.
        n = rho >= RHO(2);
    elseif iSphere == nSpheres % Capture all above top of defined atmosphere.
        n = rho < RHO(nSpheres);
    else 
        n = rho >= RHO(iSphere+1) & rho < RHO(iSphere);
    end
    
    if K(iSphere) == 0 % Isothermal layer
        h(n) = H(iSphere) - log(rho(n)/RHO(iSphere))*R*T(iSphere)/g0;
    else % Gradient layer
        TonTi = (rho(n)/RHO(iSphere)).^(-1./(1 + g0/(K(iSphere)*R)));
        Temp = TonTi*T(iSphere);
        h(n) = H(iSphere) + (Temp - T(iSphere)) / K(iSphere);
    end
end

%% Process outputs:
if dimVarOut
    h = h*u.m;
elseif convertUnits
    h = h / 0.3048;
end

end
