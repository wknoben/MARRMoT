function [out] = percolation_6(p1,p2,S,dt)
%percolation_6 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Threshold-based percolation from a store that can reach negative values
% Constraints:  f <= S/dt
% @(Inputs):    p1   - maximum percolation rate
%               p2   - storage threshold for reduced percolation [mm]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = min(S/dt,p1.*min(1,max(0,S)./p2));

end

