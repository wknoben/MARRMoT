function [out] = recharge_7(p1,fin)
%recharge_7 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Constant recharge limited by incoming flux
% Constraints:  -
% @(Inputs):    p1   - maximum recharge rate [mm/d]
%               fin  - incoming flux [mm/d]

out = min(p1,fin);

end

