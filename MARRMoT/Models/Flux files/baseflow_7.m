function [out] = baseflow_7(p1,p2,S,dt)
%baseflow_7 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear outflow from a reservoir
% Constraints:  f <= S/dt
%               S >= 0
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

out = min(S/dt,p1.*max(0,S).^p2);

end

