function [func] = evap_7(~)
%evap_7 Creates function for evaporation: evaporates based on scaled
%current water storage, limited by potential rate.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation scaled by relative storage
% Constraints:  f <= S/dt
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(S,Smax,Ep,dt) min(S./Smax.*Ep,S/dt);

end

