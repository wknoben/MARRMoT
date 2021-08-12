function [out] = evap_22(p1,p2,S,Ep,dt)
%evap_22 3-part piece-wise evaporation

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Threshold-based evaporation rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - wilting point [mm]
%               p2   - 2nd (lower) threshold [mm]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]

out = min(S/dt,max(0,min((S-p1)./(p2-p1).*Ep,Ep)));

end

