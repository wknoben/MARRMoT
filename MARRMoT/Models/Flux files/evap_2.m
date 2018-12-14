function [func] = evap_2(~)
%evap_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation at a scaled, plant-controlled rate
% Constraints:  f <= Ep
%               f <= S/dt
% @(Inputs):    p1   - plant-controlled base evaporation rate [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S,Smax,Ep,dt) min([p1*S/Smax,Ep,S/dt]);

end

