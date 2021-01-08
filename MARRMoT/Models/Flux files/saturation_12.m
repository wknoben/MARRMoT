function [out] = saturation_12(p1,p2,In)
%saturation_12 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (min-max linear variant)
% Constraints:  f >= 0
% @(Inputs):    p1   - maximum contributing fraction area [-]
%               p2   - minimum contributing fraction area [-]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = max(0,(p1-p2)./(1-p2)).*In;

end

