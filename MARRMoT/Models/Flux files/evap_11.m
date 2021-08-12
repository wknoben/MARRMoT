function [out] = evap_11(S,Smax,Ep)
%evap_11 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Evaporation quadratically related to current soil moisture
% Constraints:  f >= 0
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]

out = max(0,(2*S/Smax-(S/Smax)^2)*Ep);

end

