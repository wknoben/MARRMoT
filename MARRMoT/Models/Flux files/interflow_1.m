function [out] = interflow_1(p1,S,Smax,flux)
%interflow_1 interflow based on incoming flux size

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Interflow as a scaled fraction of an incoming flux
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               flux - incoming flux [mm/d]

out = p1*S/Smax*flux;

end

