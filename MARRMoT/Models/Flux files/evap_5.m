function [func] = evap_5(~)
%evap_5 Creates function for evaporation: evaporates based on scaled
%current water storage, for a fraction of the surface
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation from bare soil scaled by relative storage
% Constraints:  Ea <= Ep
%               Ea <= S/dt
% @(Inputs):    p1   - fraction of area that is bare soil [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S,Smax,Ep,dt) min((1-p1).*S./Smax.*Ep,S/dt);

end

