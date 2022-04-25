function [out] = abstraction_1(p1)
%abstraction_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Function for constant groundwater abstraction
% Constraints:  None, taken from a store with possible negative depth
% @(Inputs):    p1 - Abstraction rate [mm/d]

out = p1;

end

