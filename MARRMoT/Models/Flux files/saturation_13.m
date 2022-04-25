function [out] = saturation_13(p1,p2,S,In)
%saturation_13 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (normal distribution variant)
% Constraints:  -
% @(Inputs):    p1   - soil depth where 50% of catchment contributes to overland flow [mm]
%               p2   - soil depth where 16% of catchment contributes to overland flow [mm]
%               S    - current storage [mm]
%               In   - incoming flux [mm/d]

out = In.*normcdf(log10(max(0,S)./p1)./log10(p1./p2));

end

