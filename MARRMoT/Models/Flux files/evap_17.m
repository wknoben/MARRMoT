function [out] = evap_17(p1,S,Ep)
%evap_17 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Scaled evaporation from a store that allows negative values
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [mm-1]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]

out = 1./(1+exp(-1.*p1.*S)).*Ep;

end

