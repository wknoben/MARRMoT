function [out] = saturation_7(p1,p2,p3,p4,p5,S,In)
%saturation_7 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (gamma function variant)
% Constraints:  f = 0, for x-p3 < 0
%               S >= 0      prevents numerical problems with integration
% @(Inputs):    p1   - scaling parameter [-]
%               p2   - gamma function parameter [-]
%               p3   - storage threshold for flow generation [mm]
%               p4   - absolute scaling parameter [mm]
%               p5   - linear scaling parameter [-]
%               S    - current storage [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = integral(@(x)...
        1./(p1.*gamma(p2)).*(max(x-p3,0)./p1).^(p2-1).*exp(-1.*max(x-p3,0)./p1),...
        p5.*max(S,0)+p4,Inf).*In;

end

