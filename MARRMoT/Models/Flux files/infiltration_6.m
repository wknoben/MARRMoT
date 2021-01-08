function [out] = infiltration_6(p1,p2,S,Smax,fin)
%infiltration_6 Creates scaled, non-linear infiltration
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Infiltration rate non-linearly scaled by relative storage
% Constraints:  f <= fin
% @(Inputs):    p1   - base infiltration rate [mm/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               fin  - incoming flux [mm/d]
%
% WK, 07/10/2018

out = min([fin,p1.*max(0,S/Smax).^(p2).*fin]);

end

