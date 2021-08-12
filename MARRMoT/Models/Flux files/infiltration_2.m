function [out] = infiltration_2(p1,p2,S1,S1max,flux,S2,dt)
%infiltration_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  0 <= f <= S2/dt
% @(Inputs):    p1    - maximum infiltration rate [mm,/d]
%               p2    - exponential scaling parameter [-]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               flux  - reduction of infiltration rate by infiltration 
%                       demand already fulfilled elsewhere [mm/d]
%               S2    - storage available for infiltration [mm]
%               dt    - time step size [d]

out = max(min((p1.*exp(-1*p2*S1./S1max))-flux,S2/dt),0);

end

