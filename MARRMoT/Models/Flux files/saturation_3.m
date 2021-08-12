function [out] = saturation_3(S,Smax,p1,In)
%saturation_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (exponential variant)
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               p1   - linear scaling parameter [-]
%               In   - incoming flux [mm/d]

out = (1-(1/(1+exp((S/Smax + 0.5)/p1)))).*In;

end

