function [out] = baseflow_9(p1,p2,S)
%baseflow_9 Linear flow above a threshold

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Linear outflow from a reservoir if a storage threshold is exceeded
% Constraints:  -
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - storage threshold for flow generation [mm]
%               S    - current storage [mm]

out = p1.*max(0,S-p2);

end

