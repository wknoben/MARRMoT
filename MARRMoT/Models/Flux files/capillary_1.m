function [out] = capillary_1(p1,S1,S1max,S2,dt)
%capillary_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Capillary rise: based on deficit in higher reservoir
% Constraints:  f <= S2/dt
% @(Inputs):    p1   - maximum capillary rise rate  [mm/d]
%               S1   - current storage in receiving store [mm]
%               S1max- maximum storage in receiving store [mm]
%               S2   - current storage in providing store [mm]
%               dt   - time step size [d]

out = min(p1.*(1-S1/S1max),S2/dt);

end

