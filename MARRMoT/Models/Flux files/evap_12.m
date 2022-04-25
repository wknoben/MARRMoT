function [out] = evap_12(S,p1,Ep)
%evap_11 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Evaporation from deficit store, with exponential decline as
%               deficit goes below a threshold
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               p1   - wilting point [mm]
%               Ep   - potential evapotranspiration rate [mm/d]

out = min(1,exp(2*(1-S/p1)))*Ep;

end

