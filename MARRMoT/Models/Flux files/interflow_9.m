function [out] = interflow_9(S,p1,p2,p3,dt)
%interflow_9 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Non-linear interflow if storage exceeds a threshold
% Constraints:  f <= S-p2
%               S-p2 >= 0     prevents numerical issues with complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - storage threshold for flow generation [mm]
%               p3   - exponential scaling parameter [-]    
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min(max((S-p2)/dt,0),(p1.*max(S-p2,0)).^p3);

end

