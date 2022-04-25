function [out] = evap_13(p1,p2,Ep,S,dt)
%evap_13 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Exponentially scaled evaporation
% Constraints:  f <= S/dt
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               Ep   - potential evapotranspiration rate [mm/d]
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min((p1^p2)*Ep,S/dt);

end

