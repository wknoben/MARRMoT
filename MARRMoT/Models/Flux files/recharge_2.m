function [out] = recharge_2(p1,S,Smax,flux)
%recharge_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Recharge as non-linear scaling of incoming flux
% Constraints:  S >= 0
% @(Inputs):    p1   - recharge scaling non-linearity [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               flux - incoming flux [mm/d]

out = flux*((max(S,0)/Smax)^p1);

end

