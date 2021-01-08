function [out] = interception_4(p1,p2,t,tmax,In,dt)
%interception_4 Creates function for effective rainfall through interception
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Interception excess after a time-varying fraction is intercepted
% Constraints:  f >= 0
% @(Inputs):    p1   - mean throughfall fraction [-]
%               p2   - timing of maximum throughfall fraction [d]
%               t    - current time step [-]
%               tmax - duration of the seasonal cycle [d]
%               In   - incoming flux [mm/d]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = max(0,p1+(1-p1)*cos(2*pi*(t*dt-p2)/tmax))*In;

end

