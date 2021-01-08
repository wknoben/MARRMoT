function [out] = saturation_4(S,Smax,In)
%saturation_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (quadratic variant)
% Constraints:  f >= 0
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum storage [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = max(0,(1-(S/Smax).^2).*In);

end

