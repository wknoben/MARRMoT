function [func] = interception_5(~)
%interception_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Interception excess after a combined absolute amount and fraction are intercepted
% Constraints:  f >= 0
% @(Inputs):    p1   - fraction that is not throughfall [-]
%               p2   - constnat interception and evaporation [mm/d]
%               In   - incoming flux [mm/d]
%
% WK, 08/10/2018

func = @(p1,p2,In) max(p1.*In-p2,0);


end

