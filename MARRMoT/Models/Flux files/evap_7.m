function [out] = evap_7(S,Smax,Ep,dt)
%evap_7 evaporation based on scaled current water storage, limited by
%potential rate.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Evaporation scaled by relative storage
% Constraints:  f <= S/dt
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]

out = min(S./Smax.*Ep,S/dt);

end

