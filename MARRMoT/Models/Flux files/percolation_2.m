function [out] = percolation_2(p1,S,Smax,dt)
%percolation_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Percolation scaled by current relative storage
% Constraints:  f <= S/dt
% @(Inputs):    p1   - maximum percolation rate [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = min(S/dt,p1.*S/Smax);

end

