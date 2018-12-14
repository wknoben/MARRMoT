function [func] = recharge_3(~)
%recharge_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Linear recharge
% Constraints:  -
% @(Inputs):    p1   - time coefficient [d-1]
%               S    - current storage [mm]
%
% WK, 08/10/2018

func = @(p1,S) p1.*S;

end

