function [func] = baseflow_5(~)
%baseflow_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear scaled outflow from a reservoir
% Constraints:  f <= S/dt
% @(Inputs):    p1   - base outflow rate [mm/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,S,Smax,dt) min(S/dt,p1*((max(S,0)/Smax)^p2));

end

