function [out] = evap_1(S,Ep,dt)
%evap_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Evaporation at the potential rate
% Constraints:  f <= S/dt
% @(Inputs):    S    - current storage [mm]
%               Ep   - potential evaporation rate [mm/d]
%               dt   - time step size

out = min(S/dt,Ep);

end

