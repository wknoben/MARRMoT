function [func] = evap_9(~)
%evap_9 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation from bare soil scaled by relative storage and 
%               by relative water availability across all stores
%               all stores
% Constraints:  f <= S/dt
%               f >= 0
% @(Inputs):    S1   - current storage in store 1 [mm]
%               S2   - current storage in store 2 [mm]
%               p1   - fraction vegetated area [-]
%               p2   - wilting point [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(S1,S2,p1,Smax,Ep,dt) max(min(S1/(S1+S2)*(1-p1)*S1/(Smax-S2)*Ep,S1/dt),0);

end

