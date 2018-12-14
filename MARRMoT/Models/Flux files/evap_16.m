function [func] = evap_16(~)
%evap_16 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Scaled evaporation if another store is below a threshold
% Constraints:  f <= S1/dt
% @(Inputs):    p1    - linear scaling parameter [-]
%               S1    - current storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2min - threshold S2 storage for evaporation occurence [mm]
%               Ep    - potential evapotranspiration rate [mm/d]
%               dt    - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S1,S2,S2min,Ep,dt) min((p1.*Ep).*smoothThreshold_storage_logistic(S2,S2min),S1/dt);

end

