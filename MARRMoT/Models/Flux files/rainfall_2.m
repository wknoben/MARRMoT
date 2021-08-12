function [out] = rainfall_2(In,T,p1,p2)
%rainfall_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Rainfall based on a temperature threshold interval
% Constraints:  -
% @(Inputs):    In   - incoming precipitation flux [mm/d]
%               T    - current temperature [oC]
%               p1   - midpoint of the combined rain/snow interval [oC]
%               p2   - length of the mixed snow/rain interval [oC]

out = min(In,max(0,In.*(T-(p1-0.5*p2))/p2));

end

