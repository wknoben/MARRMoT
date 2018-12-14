function [func] = routing_1(~)
%routing_1 Creates function for non-linear routing
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Threshold-based non-linear routing
% Constraints:  f <= S/dt
%               S >= 0      prevents complex numbers
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               p3   - fractional release parameter [-]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(p1,p2,p3,S,dt) min([S/dt,p1.*(max(S,0).^p2),p3.*S/dt]);

end

