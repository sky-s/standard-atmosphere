function h = densityalt(rho,varargin)
% DENSITYALT  Returns altitude corresponding to the given array of air densities
%   in the standard or non-standard atmosphere.
%
%   H = DENSITYALT(RHO) returns altitude, h, as a function of air density, rho.
% 
%   The input RHO can be followed by parameter/value pairs for further control
%   of DENSITYALT. Possible parameters are:
%     inputUnits     - String for units of input RHO, either kg/m³ or slug/ft³.
%                      [{SI}|kg/m3|kg/m^3  |  US|slug/ft3|slug/ft^3]
%     outputUnits    - String for units of output H, either meters or feet.
%                      [{SI}|m|meters  |  US|ft|feet]
%     atmosphereFunc - String determining atmosphere function to be used.
%                      [{atmos}|tropos|atmosisa|atmoscoesa|atmosnonstd]
%     atmosphereArgs - Cell array of additional arguments to pass to 
%                      atmosphereFunction after the density input (e.g. for non-
%                      standard atmospheres).
%     method         - Method used for either searching for or interpolating a
%                      solution (the equations that define the standard 
%                      atmosphere cannot be inverted in terms of density).
%                      Search: [fzero | bisection] (fzero only for scalar case)
%                      Interpolate: any method accepted by interp1.
%                      Default method is 'pchip' for interpolation.
%     options        - Options used for fzero or bisection methods.
%     hMin           - For search: lower search interval bound in meters. 
%                      For interpolation: start of generated interpolation grid.
%                      Default hMin = 0.
%     hMax           - For search: upper search interval bound in meters. 
%                      For interpolation: end of generated interpolation grid.
%                      Default hMax = 86000.
%     spacing        - Spacing of of interpolation grid in meters. 
%                      Default spacing = 50.
% 
%   If the input RHO is a DimVar, inputUnits and outputUnits will be ignored and
%   the output will be a DimVar.
%
%   DENSITYALT is valid for the entire standard atmosphere up through the
%   mesopause (86 km height). It assumes that all that is known is air
%   density. If pressure or temperature and density are known, there exist
%   more straightforward methods for calculating density altitude:
%       P = rho*R*T; h = h0 * (1 - P^0.190284), where P is in atmospheres and 
%           h0 = 145366.45 ft or 44307.694 m.
%       (http://www.srh.noaa.gov/images/epz/wxcalc/pressureAltitude.pdf)
% 
%   Example: Plot altitude as a function of air density for a cold (-15°C),
%   standard, and hot (+15°C) day (leverages vectorization of bisection).
%       rho = 1.225:-0.01:0.025; % kg/m^3
%       tempOffset = [-15 0 15]; % °C
%       [rho, tempOffset] = meshgrid(rho,tempOffset);
%       h = densityalt(rho,'method','bisection','outputUnits','ft',...
%               'atmosphereArgs',{tempOffset},'hMin',-5000);
%       plot(rho',h'/1000);
%       xlabel('Air density (kg/m^3)'); ylabel('Altitude (kft)'); grid on
%       legend('Cold','Standard','Hot')
% 
%   See also ATMOSISA, ATMOSNONSTD, ATMOSCOESA, INTERP1, FZERO, 
%     BISECTION - http://www.mathworks.com/matlabcentral/fileexchange/28150,
%     ATMOS     - http://www.mathworks.com/matlabcentral/fileexchange/28135,
%     U, UNITS  - http://www.mathworks.com/matlabcentral/fileexchange/38977.
%
%   H = DENSITYALT(RHO,Param1,Val1,Param2,Val2,...)

%   Copyright 2012-2013, 2015 Sky Sartorius
%   Author contact: mathworks.com/matlabcentral/fileexchange/authors/101715

% The 'pchip' method is smooth, preserves the shape, and has very tiny errors
% for most of the domain. However, the errors around the discontinuous
% transisions between atmospheres (e.g. between troposphere and tropopause) can
% be on the order of 1/100th of the spacing of given H. For spacing of 50 m,
% errors are never more than 813 mm. The absolute errors at these
% discontinuities are an order of magnitude lower for the 'linear' method but
% the error is present across the entire domain.

%% Parse inputs:
p = inputParser;
fName = 'densityalt';
p.FunctionName = fName;

p.addRequired('rho',@(x) validateattributes(x,{'numeric','DimVar'},...
    {'positive'},'','density, rho,'));

% Ignored if input is DimVar:
p.addParameter('inputUnits','SI');
p.addParameter('outputUnits','SI');

p.addParameter('atmosphereFunction','atmos');
p.addParameter('atmosphereArgs',{},@iscell);

p.addParameter('method','pchip'); 
% Possibilities are exact matches of 'fzero' or 'bisection', or an interp1
% method. Validation by interp1.

p.addParameter('hMin',0,@(x) validateattributes(x,{'numeric'},...
    {'scalar','finite','real'},'','altitude minimum h (m)'));
p.addParameter('hMax',86000,@(x) validateattributes(x,{'numeric'},...
    {'scalar','finite','real'},'','altitude maximum h (m)'));

% Used only with fzero/bisection search method:
p.addParameter('options',{})

% Used only with an interp method:
p.addParameter('spacing',50,@(x) validateattributes(x,{'numeric'},...
    {'positive','scalar','finite','real'},'','altitude grid spacing (m)'));


parse(p,rho,varargin{:});
i = p.Results;
rho = i.rho;
if ~iscell(i.options)
    i.options = {i.options};
end

% Validate inputs strings:
i.inputUnits = validatestring(i.inputUnits,...
    {'SI','kg/m3','kg/m^3','US','slug/ft3','slug/ft^3'},fName,'inputUnits');
i.outputUnits = validatestring(i.outputUnits,...
    {'SI','m','meters','US','ft','feet'},fName,'outputUnits');

i.atmosphereFuntion = validatestring(i.atmosphereFunction,...
    {'atmos','tropos','atmosisa','atmoscoesa','atmosnonstd'},...
    fName,'atmosphereFunction');


%% Process input density and return density in units of kg/m^3:
dimVarOut = false; % Flag to convert output to a DimVar.
if isa(rho,'DimVar')
    rho = rho/(u.kg/u.m^3);
    dimVarOut = true; 
elseif any(strcmpi(i.inputUnits,{'US','slug/ft3','slug/ft^3'}))
    % Convert from imperial units.
    rho = rho * 515.3788183931961;
end
% Otherwise, input rho is already in kg/m^3.

%% Build atmosphere model:
% Build an atmosphere function that takes only altitude (in meters) as input and
% returns density (in kg/m^3) as the first and only output.
function rho = myAtmo(h)
    switch i.atmosphereFunction
        case {'atmos' 'tropos'}
            rho = feval(i.atmosphereFunction,h,i.atmosphereArgs{:});
        otherwise
            % Other functions return density as the fourth output.
            [~,~,~,rho] = feval(i.atmosphereFunction,h,i.atmosphereArgs{:});
    end
end

%% Find the altitude:

% Switch between search and interpolation methods:
switch i.method
    case 'fzero'
        h = fzero(@(h) rho-myAtmo(h),[i.hMin,i.hMax],i.options{:});
        
    case 'bisection'
        h = bisection(@myAtmo,i.hMin,i.hMax,rho,i.options{:});
        
    otherwise % interp1
        H = (i.hMin:i.spacing:i.hMax)';
        RHO = myAtmo(H);
        h = interp1(RHO,H,rho,i.method);
end


%% Process output altitude from meters into desired units:
if dimVarOut
    h = h*u.m;
elseif any(strcmpi(i.outputUnits,{'US','ft','feet'})) 
    % Convert to imperial units if necessary.
    h = h / 0.3048;
end
% Otherwise, output h is already in meters.

end