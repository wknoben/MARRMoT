function [out] = recharge_1(p1,S,Smax,flux)
%recharge_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Recharge as scaled fraction of incoming flux
% Constraints:  -
% @(Inputs):    p1   - fraction of flux that is recharge [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               flux - incoming flux [mm/d]

out = p1*S/Smax*flux;

end

