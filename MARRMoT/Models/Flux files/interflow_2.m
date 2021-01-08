function [out] = interflow_2(p1,S,p2,dt)
%interflow_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear interflow
% Constraints:  f <= S
%               S >= 0  - this avoids numerical issues with complex numbers
% @(Inputs):    p1   - time delay [d-1]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = min(p1*max(S,0)^(1+p2),max(S/dt,0));

end

