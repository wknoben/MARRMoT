function [out] = interception_3(p1,In)
%interception_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Interception excess after a fraction is intercepted
% Constraints:  -
% @(Inputs):    p1   - fraction throughfall [-]
%               In   - incoming flux [mm/d]
%
% WK, 08/10/2018

out = p1.*In;

end

