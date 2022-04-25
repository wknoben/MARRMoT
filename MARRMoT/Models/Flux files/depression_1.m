function [out] = depression_1(p1,p2,S,Smax,flux,dt)
%depression_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Exponential inflow to surface depression store
% Constraints:  f <= (Smax-S)/dt
%               S <= Smax
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]

out = min(p1.*exp(-1.*p2.*S./max(Smax-S,0)).*flux,max((Smax-S)/dt,0));

end

