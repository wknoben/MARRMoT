function [out] = interflow_7(S,Smax,p1,p2,p3,dt)
%interflow_7 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Non-linear interflow if storage exceeds a threshold
% Constraints:  f <= (S-p1*Smax)/dt
%               S-p1*Smax >= 0      prevents numerical issues with complex
%                                   numbers
% @(Inputs):    p1   - storage threshold as fraction of Smax [-]
%               p2   - time coefficient [d]
%               p3   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]

out = min(max(0,(S-p1.*Smax)/dt),(max(0,S-p1.*Smax)./p2).^(1/p3));

end

