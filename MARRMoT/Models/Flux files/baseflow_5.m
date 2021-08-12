function [out] = baseflow_5(p1,p2,S,Smax,dt)
%baseflow_5 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Non-linear scaled outflow from a reservoir
% Constraints:  f <= S/dt
% @(Inputs):    p1   - base outflow rate [mm/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               dt   - time step size [d]

out = min(S/dt,p1*((max(S,0)/Smax)^p2));

end

