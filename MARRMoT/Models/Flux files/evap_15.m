function [out] = evap_15(Ep,S1,S1max,S2,S2min,dt)
%evap_15 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Scaled evaporation if another store is below a threshold
% Constraints:  f <= S1/dt
% @(Inputs):    Ep    - potential evapotranspiration rate [mm/d]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]
%               dt    - time step size [d]
%
% WK, 05/10/2018

out = min((S1/S1max*Ep).*smoothThreshold_storage_logistic(S2,S2min,S1/dt));

end

