function [func] = baseflow_6(~)
%baseflow_6 Creates function for quadratic baseflow
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Quadratic outflow from a reservoir if a storage threshold is exceeded
% Constraints:  f <= S/dt
% @(Inputs):    p1   - linear scaling parameter [mm-1 d-1]
%               p2   - threshold that must be exceeded for flow to occur [mm]
%               S    - current storage [mm]
%
% WK, 05/10/2018

func = @(p1,p2,S,dt) min(S/dt,p1.*S.^2).*(1-smoothThreshold_storage_logistic(S,p2));

end

