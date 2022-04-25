function [out] = saturation_12(p1,p2,In)
%saturation_12 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (min-max linear variant)
% Constraints:  f >= 0
% @(Inputs):    p1   - maximum contributing fraction area [-]
%               p2   - minimum contributing fraction area [-]
%               In   - incoming flux [mm/d]

out = max(0,(p1-p2)./(1-p2)).*In;

end

