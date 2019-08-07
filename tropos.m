function [rho,a,temp,press,kvisc,sigma]=tropos(h_in,tOffset)
% TROPOS  Stripped-down version of atmos, applicable only to the troposphere
% (covers the vast majority of atmospheric flight), for when computation speed
% is a priority.
% 
%   [rho,a,T,P,nu,sigma] = TROPOS(h)
%   [rho,a,T,P,nu,sigma] = TROPOS(h,dT)
%   
%   See also ATMOS.


if nargin < 2
    tOffset = 0;
end
if nargin < 1
    h_in = 0;
end

dimVarOut = false;
if isa(h_in,'DimVar')
    h_in = h_in/u.m;
    dimVarOut = true;
end
if isa(tOffset,'DimVar')
    tOffset = tOffset/u.K;
    % It is allowed to mix DimVar h_in and double tOffset (or reverse). 
end

% h_in(h_in>11000 | h_in<0) = NaN;

TonTi=1-2.255769564462953e-005*h_in;
press=101325*TonTi.^5.255879734954165;
temp = TonTi*288.15 + tOffset;
rho = press./temp/287.0528742470439;
sigma = rho/1.225;

a = sqrt(401.874018 * temp);
kvisc = (1.458e-6 * temp.^1.5 ./ (temp + 110.4)) ./ rho;

if dimVarOut
    rho = rho*u.kg/(u.m^3);
    a = a*u.m/u.s;
    temp = temp*u.K;
    press = press*u.Pa;
    kvisc = kvisc*u.m^2/u.s;
end