function [out] = evap_2(p1,S,Smax,Ep,dt)
%evap_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Evaporation at a scaled, plant-controlled rate
% Constraints:  f <= Ep
%               f <= S/dt
% @(Inputs):    p1   - plant-controlled base evaporation rate [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]

out = min([p1*S/Smax,Ep,S/dt]);

end

