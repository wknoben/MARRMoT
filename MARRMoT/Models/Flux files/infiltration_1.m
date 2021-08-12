function [out] = infiltration_1(p1,p2,S,Smax,fin)
%infiltration_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  f <= fin
% @(Inputs):    p1   - maximum infiltration rate [mm,/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]

out = min(p1.*exp((-1*p2*S)./Smax),fin);

end

