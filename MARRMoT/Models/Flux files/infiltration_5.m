function [out] = infiltration_5(p1,p2,S1,S1max,S2,S2max)
%infiltration_5 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Maximum infiltration rate non-linearly based on relative deficit and storage
% Constraints:  S2 >= 0     - prevents complex numbers
%               f <= 10^9   - prevents numerical issues with Inf outcomes
% @(Inputs):    p1    - base infiltration rate [mm/d]
%               p2    - exponential scaling parameter [-]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]

out = max(0,min(10^9,p1.*(1-S1./S1max).*max(0,S2./S2max).^(-1.*p2)));

end

