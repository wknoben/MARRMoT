function [out] = snowfall_2(In,T,p1,p2)
%snowfall_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Snowfall based on a temperature threshold interval
% Constraints:  -
% @(Inputs):    p1   - midpoint of the combined rain/snow interval [oC]
%               p2   - length of the mixed snow/rain interval [oC]
%               T    - current temperature [oC]
%               In   - incoming precipitation flux [mm/d]

out = min(In,max(0,In.*(p1+0.5*p2-T)/p2));

end

