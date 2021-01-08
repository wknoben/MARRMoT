function [out] = percolation_3(S,Smax)
%percolation_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear percolation (empirical)
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%
% WK, 08/10/2018

out = Smax^(-4)./(4)*(4/9)^(-4)*S^5;

end

