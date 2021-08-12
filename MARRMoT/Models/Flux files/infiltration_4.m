function [out] = infiltration_4(fluxIn,p1)
%infiltration_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Constant infiltration rate 
% Constraints:  f <= fin
% @(Inputs):    p1   - Infiltration rate [mm/d]
%               fin  - incoming flux [mm/d]

out = min(fluxIn,p1);

end

