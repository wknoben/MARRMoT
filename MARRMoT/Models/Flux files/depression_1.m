function [func] = depression_1(~)
%depression_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Exponential inflow to surface depression store
% Constraints:  f <= (Smax-S)/dt
%               S <= Smax
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,S,Smax,flux,dt) min(p1.*exp(-1.*p2.*S./max(Smax-S,0)).*flux,max((Smax-S)/dt,0));

end

