function [out] = baseflow_1(p1,S)
% baseflow_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Outflow from a linear reservoir
% Constraints:  -
% @(Inputs):    p1   - time scale parameter [d-1]
%               S    - current storage [mm]

out = p1.*S;

end