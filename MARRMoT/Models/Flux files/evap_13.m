function [func] = evap_13(~)
%evap_13 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Exponentially scaled evaporation
% Constraints:  f <= S/dt
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               Ep   - potential evapotranspiration rate [mm/d]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,Ep,S,dt) min((p1^p2)*Ep,S/dt);

end

