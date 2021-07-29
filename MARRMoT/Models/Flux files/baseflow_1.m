function [out] = baseflow_1(p1,S)
%baseflow_1 
%
% Copyright (C) 2021 L. Trotter
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Flux function
% ------------------
% Description:  Outflow from a linear reservoir
% Constraints:  -
% @(Inputs):    p1   - time scale parameter [d-1]
%               S    - current storage [mm]

out = p1.*S;

end