function [func] = exchange_1(~)
%exchange_1 Creates function for two-way channel exchange: linear and exponential.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Water exchange between aquifer and channel
% Constraints:  f <= fIn
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - linear scaling parameter [-]
%               p3   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               fmax - maximum flux size [mm/d]
%               dt   - time step size [d]
%
% WK, 07/10/2018

func = @(p1,p2,p3,S,fmax,dt) max((p1*abs(S/dt) + p2*(1-exp(-1*p3*abs(S/dt)))).*sign(S),-1*fmax);

end

