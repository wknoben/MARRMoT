function [out] = percolation_2(p1,S,Smax,dt)
%percolation_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Percolation scaled by current relative storage
% Constraints:  f <= S/dt
% @(Inputs):    p1   - maximum percolation rate [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]

out = min(S/dt,p1.*S/Smax);

end

