function [out] = routing_1(p1,p2,p3,S,dt)
%routing_1 non-linear routing

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Threshold-based non-linear routing
% Constraints:  f <= S/dt
%               S >= 0      prevents complex numbers
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               p3   - fractional release parameter [-]
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min([S/dt,p1.*(max(S,0).^p2),p3.*S/dt]);

end

