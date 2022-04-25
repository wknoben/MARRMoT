function [out] = baseflow_4(p1,p2,S)
%baseflow_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Exponential outflow from deficit store
% Constraints:  - 
% @(Inputs):    p1   - base outflow rate [mm/d]
%               p2   - exponent parameter [mm-1]
%               S    - current storage [mm]

out = p1*exp(-1*p2*S);

end

