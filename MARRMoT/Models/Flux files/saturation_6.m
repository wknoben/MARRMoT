function [out] = saturation_6(p1,S,Smax,In)
%saturation_6 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (linear variant)
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               In   - incoming flux [mm/d]

out = p1*S/Smax*In;

end

