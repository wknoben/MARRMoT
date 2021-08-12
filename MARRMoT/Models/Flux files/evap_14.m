function [out] = evap_14(p1,p2,Ep,S1,S2,S2min,dt)
%evap_14 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Exponentially scaled evaporation that only activates if 
%               another store goes below a certain threshold
% Constraints:  f <= S1/dt
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               Ep   - potential evapotranspiration rate [mm/d]
%               S1   - current storage in S1 [mm]
%               S2   - current storage in S2 [mm]
%               S2min- threshold for evaporation deactivation [mm]
%               dt   - time step size [d]

out = min((p1^p2)*Ep,S1/dt).*smoothThreshold_storage_logistic(S2,S2min);

end

