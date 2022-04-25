function [out] = recharge_6(p1,p2,S,dt)
%recharge_6 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Recharge to fulfil evaporation demand if the receiving 
%               store is below a threshold
% Constraints:  f <= S/dt
%               S >= 0      prevents complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - non-linear scaling [mm]
%               S    - current storage [mm]
%               dt   - time step size [d]

out = min(max(S/dt,0),p1.*max(S,0).^p2);

end

