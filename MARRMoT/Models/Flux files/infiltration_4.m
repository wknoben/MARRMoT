function [out] = infiltration_4(fluxIn,p1)
%infiltration_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Constant infiltration rate 
% Constraints:  f <= fin
% @(Inputs):    p1   - Infiltration rate [mm/d]
%               fin  - incoming flux [mm/d]
%
% WK, 07/10/2018

out = min(fluxIn,p1);

end

