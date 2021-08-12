function [out] = interception_4(p1,p2,t,tmax,In,dt)
%interception_4 Creates function for effective rainfall through interception

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Interception excess after a time-varying fraction is intercepted
% Constraints:  f >= 0
% @(Inputs):    p1   - mean throughfall fraction [-]
%               p2   - timing of maximum throughfall fraction [d]
%               t    - current time step [-]
%               tmax - duration of the seasonal cycle [d]
%               In   - incoming flux [mm/d]
%               dt   - time step size [d]

out = max(0,p1+(1-p1)*cos(2*pi*(t*dt-p2)/tmax))*In;

end

