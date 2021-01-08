function [out] = saturation_5(S,p1,p2,In)
%saturation_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Deficit store: exponential saturation excess based on current 
%               storage and a threshold parameter
% Constraints:  S >= 0      prevents numerical issues with complex numbers
% @(Inputs):    p1   - deficit threshold above which no flow occurs [mm]
%               p2   - exponential scaling parameter [-]
%               S    - current deficit [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = (1-min(1,(max(S,0)/p1).^p2)).*In;

end

