function [out] = capillary_3(p1,p2,S1,S2,dt)
%capillary_3 capillary rise: linear relation

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Capillary rise scaled by receiving store's deficit up to a storage threshold
% Constraints:  f <= S2
% @(Inputs):    p1   - base capillary rise rate [mm/d]
%               p2   - threshold above which no capillary flow occurs [mm]
%               S1   - current storage in receiving store [mm]
%               S2   - current storage in supplying store [mm]
%               dt   - time step size [d]

out = min(S2/dt,p1*(1-S1/p2)*smoothThreshold_storage_logistic(S1,p2));

end

