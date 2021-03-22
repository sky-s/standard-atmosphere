function h = pressurealt(P,units)
% PRESSUREALT  Approximate altitude in troposphere for a given pressure.
% 
%   h = PRESSUREALT(P,units) returns altitude given provided array of pressures.
%   
%   units is a string specifying expected input and returned output units,
%   either 'SI' (default) or 'US'. This is ignored if the provided input is a
%   DimVar, in which case all outputs are also DimVars.
%                 Description:         SI:           US:
%                 --------------       -----         -----
%       h         Altitude             m             ft
%       P         Pressure             Pa            lbf/ft^2      
%   
%   See also atmos.

% References:
%   https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf
%   http://en.wikipedia.org/wiki/Pressure_altitude

if isa(P,'DimVar')
    h = u.ft * 145366.45 * (1 - (P/u.atm).^0.190284);
elseif nargin < 2 || strcmpi(units,"SI")
    h = 44307.69396 * (1 - (P/101325).^0.190284);
elseif strcmpi(units,"US")
    h = 145366.45 * (1 - (P/2116.21662367394).^0.190284);
else
    error('Unknown units.')
end
