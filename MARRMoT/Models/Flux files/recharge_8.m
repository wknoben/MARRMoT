function [out] = recharge_8(p1,S,Smax,p2,dt)
%recharge_8 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Recharge as non-linear scaling of incoming flux
% Constraints:  f <= S/dt
%               S >= 0
% @(Inputs):    p1   - recharge scaling non-linearity [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               p2   - maximum flux rate [mm/d]
%               dt   - time step size [d]

out = min(p2*((max(S,0)/Smax)^p1),max(S/dt,0));

end

