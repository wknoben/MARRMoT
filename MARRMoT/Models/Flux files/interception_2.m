function [func] = interception_2(~)
%interception_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Interception excess after a constant amount is intercepted
% Constraints:  f >= 0
% @(Inputs):    In   - incoming flux [mm/d]
%               p1   - interception and evaporation capacity [mm/d]
%
% WK, 07/10/2018

func = @(In,p1) max(In-p1,0);

end

