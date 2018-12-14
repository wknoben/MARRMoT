function [func] = recharge_8(~)
%recharge_8 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Recharge as non-linear scaling of incoming flux
% Constraints:  f <= S/dt
%               S >= 0
% @(Inputs):    p1   - recharge scaling non-linearity [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               p2   - maximum flux rate [mm/d]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(p1,S,Smax,p2,dt) min(p2*((max(S,0)/Smax)^p1),max(S/dt,0));

end

