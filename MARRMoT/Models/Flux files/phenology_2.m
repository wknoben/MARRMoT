function [out] = phenology_2(p1,p2,p3,t,tmax,dt)
%phenology_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Phenology-based maximum interception capacity (returns store size [mm])
% Constraints:  Implicit assumption: 0 <= p2 <= 1
% @(Inputs):    p1   - mean interception capacity [mm]
%               p2   - seasonal change as fraction of the mean [-]
%               p3   - time of maximum store size [d]
%               t    - current time step [-]
%               tmax - seasonal length [d]
%               dt   - time step size [d]

out = p1*(1+p2*sin(2*pi*(t*dt-p3)/tmax));

end

