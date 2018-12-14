function [func] = evap_4(~)
%evap_4 Creates function for evaporation: evaporates based on scaled
%current water storage, a walting point, a constraining factor and limited 
%by potential rate.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Constrained, scaled evaporation if storage is above a wilting point
% Constraints:  f <= S/dt
% @(Inputs):    Ep   - potential evapotranspiration rate [mm/d]
%               p1   - scaling parameter [-]
%               p2   - wilting point as fraction of Smax [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(Ep,p1,S,p2,Smax,dt) min(Ep.*max(0,p1*(S-p2.*Smax)./(Smax-p2.*Smax)),S/dt);

end

