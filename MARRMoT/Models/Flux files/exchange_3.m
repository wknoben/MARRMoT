function [out] = exchange_3(p1,S,p2)
%exchange_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Water exchange with infinite size store based on threshold
% Constraints:  -
% @(Inputs):    p1   - base leakage time delay [d-1]
%               p2   - threshold for flow reversal [mm]
%               S    - current storage [mm]

out = p1*(S-p2);

end

