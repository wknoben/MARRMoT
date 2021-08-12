function [out] = evap_18(p1,p2,p3,S,Ep)
%evap_18 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Exponentially declining evaporation from deficit store
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - linear scaling parameter [-]
%               p3   - storage scaling parameter [mm]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]

out = p1.*exp(-1.*p2.*S./p3).*Ep;

end

