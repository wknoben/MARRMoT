function [func] = evap_22(~)
%evap_22 Creates 3-part piece-wise evaporation
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Threshold-based evaporation rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - wilting point [mm]
%               p2   - 2nd (lower) threshold [mm]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,S,Ep,dt) min(S/dt,max(0,min((S-p1)./(p2-p1).*Ep,Ep)));


end

