function [func] = saturation_10(~)
%saturation_10 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (min-max exponential variant)
% Constraints:  -
% @(Inputs):    p1   - maximum contributing fraction area [-]
%               p2   - minimum contributing fraction area [-]
%               p3   - exponentia scaling parameter [-]
%               S    - current storage [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

func = @(p1,p2,p3,S,In) min(p1,p2+p2.*exp(p3.*S)).*In;


end

