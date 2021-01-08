function [out] = recharge_4(p1,S,dt)
%recharge_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Constant recharge
% Constraints:  f <= S/dt
% @(Inputs):    p1   - time coefficient [d-1]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = min(p1,S/dt);

end

