function [out] = baseflow_2(S,p1,p2,dt)
%baseflow_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear outflow from a reservoir
% Constraints:  f <= S/dt
%               S >= 0       prevents numerical issues with complex numbers
% @(Inputs):    S    - current storage [mm]
%               p1   - time coefficient [d]
%               p2   - exponential scaling parameter [-]
%               dt   - time step size [d]
%
% WK, 05/10/2018

out = min((1./p1*max(S,0)).^(1./p2),max(S,0)/dt);

end

