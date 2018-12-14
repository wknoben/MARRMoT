function [func] = evap_14(~)
%evap_14 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
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
%
% WK, 05/10/2018

func = @(p1,p2,Ep,S1,S2,S2min,dt) min((p1^p2)*Ep,S1/dt).*smoothThreshold_storage_logistic(S2,S2min);

end

