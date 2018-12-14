function [func] = percolation_1(~)
%percolation_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Percolation at a constant rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - base percolation rate [mm/d]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(p1,S,dt) min(p1,S/dt);

end

