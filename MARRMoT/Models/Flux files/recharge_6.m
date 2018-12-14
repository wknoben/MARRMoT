function [func] = recharge_6(~)
%recharge_6 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Recharge to fulfil evaporation demand if the receiving 
%               store is below a threshold
% Constraints:  f <= S/dt
%               S >= 0      prevents complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - non-linear scaling [mm]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(p1,p2,S,dt) min(max(S/dt,0),p1.*max(S,0).^p2);

end

