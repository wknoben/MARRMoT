function [out] = percolation_1(p1,S,dt)
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

out = min(p1,S/dt);

end

