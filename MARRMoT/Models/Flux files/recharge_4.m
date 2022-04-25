function [out] = recharge_4(p1,S,dt)
%recharge_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Constant recharge
% Constraints:  f <= S/dt
% @(Inputs):    p1   - time coefficient [d-1]
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min(p1,S/dt);

end

