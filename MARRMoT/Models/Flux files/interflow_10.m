function [out] = interflow_10(S,p1,p2,p3)
%interflow_10 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Scaled linear interflow if storage exceeds a threshold
% Constraints:  f = 0, for S < p2
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - threshold for flow generation [mm]
%               p3   - linear scaling parameter [-]
%               S    - current storage [mm]

out = p1*max(0,S-p2)/(p3);

end

