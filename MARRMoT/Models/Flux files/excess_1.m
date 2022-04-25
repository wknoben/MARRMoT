function [out] = excess_1(So,Smax,dt)
%excess_1

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Storage excess when store size changes (returns flux [mm/d])
% Constraints:  f >= 0
% @(Inputs):    So   - 'old' storage [mm]
%               Smax - 'new' maximum storage [mm]
%               dt   - time step size [d]

out = max((So-Smax)/dt,0);

end

