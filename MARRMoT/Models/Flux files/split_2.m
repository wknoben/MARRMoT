function [func] = split_2(~)
%split_2
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Split flow (returns flux [mm/d]), counterpart to split_1
% Constraints:  -
% @(Inputs):    p1   - fraction of flux to be diverted [-]
%               In   - incoming flux [mm/d]

func = @(p1,In) (1-p1).*In;

end