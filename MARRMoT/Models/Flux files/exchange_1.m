function [out] = exchange_1(p1,p2,p3,S,fmax,dt)
%exchange_1 two-way channel exchange: linear and exponential.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Water exchange between aquifer and channel
% Constraints:  f <= fIn
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - linear scaling parameter [-]
%               p3   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               fmax - maximum flux size [mm/d]
%               dt   - time step size [d]

out = max((p1*abs(S/dt) + p2*(1-exp(-1*p3*abs(S/dt)))).*sign(S),-1*fmax);

end

