function [func] = split_1(~)
%split_1 Creates function for flow splitting
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Split flow (returns flux [mm/d])
% Constraints:  -
% @(Inputs):    p1   - fraction of flux to be diverted [-]
%               In   - incoming flux [mm/d]
%
% WK, 08/10/2018

func = @(p1,In) p1.*In;

end

