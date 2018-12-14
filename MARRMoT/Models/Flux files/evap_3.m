function [func] = evap_3(~)
%evap_3 Creates function for evaporation: 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation based on scaled current water storage and wilting point
% Constraints:  f <= Ep
%               f <= S/dt
% @(Inputs):    p1   - wilting point as fraction of Smax [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S,Smax,Ep,dt) min([S/(p1*Smax)*Ep,Ep,S/dt]);

end

