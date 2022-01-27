function [out] = percolation_3(S,Smax)
%percolation_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Non-linear percolation (empirical)
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]

out = Smax^(-4)./(4)*(4/9)^(4)*S^5;

end

