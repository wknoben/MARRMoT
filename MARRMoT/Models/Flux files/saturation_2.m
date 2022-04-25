function [out] = saturation_2(S,Smax,p1,In)
%saturation_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess from a store with different degrees of saturation
% Constraints:  1-S/Smax >= 0       prevents numerical issues with complex
%                                   numbers
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               p1   - non-linear scaling parameter [-]
%               In   - incoming flux [mm/d]

% NOTE: When stores are very slightly below or over their maximum, the
% exponent can push this function into regions where no feasible solutions
% exist. The min(max()) combination prevents this from happening. 

out = (1- min(1,max(0,(1-S./Smax))).^p1) .*In;

end

