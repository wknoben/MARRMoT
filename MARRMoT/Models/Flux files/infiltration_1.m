function [out] = infiltration_1(p1,p2,S,Smax,fin)
%infiltration_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  f <= fin
% @(Inputs):    p1   - maximum infiltration rate [mm,/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]
%
% WK, 07/10/2018

out = min(p1.*exp((-1*p2*S)./Smax),fin);

end

