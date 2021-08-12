function [out] = interflow_3(p1,p2,S,dt)
%interflow_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Non-linear interflow (variant)
% Constraints:  f <= S
%               S >= 0  - this avoids numerical issues with complex numbers
% @(Inputs):    p1   - time delay [d-1]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min(p1*max(S,0)^(p2),max(S/dt,0));

end

