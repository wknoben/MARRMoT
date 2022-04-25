function [out] = baseflow_6(p1,p2,S,dt)
%baseflow_6 Creates function for quadratic baseflow

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Quadratic outflow from a reservoir if a storage threshold is exceeded
% Constraints:  f <= S/dt
% @(Inputs):    p1   - linear scaling parameter [mm-1 d-1]
%               p2   - threshold that must be exceeded for flow to occur [mm]
%               S    - current storage [mm]

out = min(S/dt,p1.*S.^2).*(1-smoothThreshold_storage_logistic(S,p2));

end

