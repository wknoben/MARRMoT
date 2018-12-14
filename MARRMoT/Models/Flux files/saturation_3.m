function [func] = saturation_3(~)
%saturation_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (exponential variant)
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               p1   - linear scaling parameter [-]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

func = @(S,Smax,p1,In) (1-(1/(1+exp((S/Smax + 0.5)/p1)))).*In;

end

