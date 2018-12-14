function [func] = evap_19(~)
%evap_19 Creates scaled, non-linear evaporation
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear scaled evaporation
% Constraints:  f <= Ep
%               f <= S/dt
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,S,Smax,Ep,dt) min([S/dt,Ep,p1.*max(0,S/Smax).^(p2).*Ep]);

end

