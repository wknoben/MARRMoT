function [out] = saturation_6(p1,S,Smax,In)
%saturation_6 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (linear variant)
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = p1*S/Smax*In;

end

