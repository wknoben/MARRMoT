function [func] = interception_3(~)
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

func = @(p1,In) p1.*In;

end

