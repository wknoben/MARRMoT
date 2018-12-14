function [func] = capillary_2(~)
%capillary_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Capillary rise at constant rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - base capillary rise rate [mm/d]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S,dt) min(p1,S/dt);

end

