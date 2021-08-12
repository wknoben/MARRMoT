function [out] = evap_4(Ep,p1,S,p2,Smax,dt)
%evap_4 evaporation based on scaled current water storage, a walting point,
%a constraining factor and limited by potential rate.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Constrained, scaled evaporation if storage is above a wilting point
% Constraints:  f <= S/dt
% @(Inputs):    Ep   - potential evapotranspiration rate [mm/d]
%               p1   - scaling parameter [-]
%               p2   - wilting point as fraction of Smax [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]

out = min(Ep.*max(0,p1*(S-p2.*Smax)./(Smax-p2.*Smax)),S/dt);

end

