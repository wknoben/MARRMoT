function [out] = evap_21(p1,p2,S,Ep,dt)
%evap_21 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Threshold-based evaporation with constant minimum rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - wilting point (1st threshold) [mm]
%               p2   - 2nd threshold as fraction of wilting point [-]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]

out = min(max(p2,min(S./p1,1)).*Ep,S./dt);

end

