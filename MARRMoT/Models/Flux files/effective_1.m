function [out] = effective_1(In1,In2)
%effective_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  General effective flow (returns flux [mm/d])
% Constraints:  In1 > In2
% @(Inputs):    In1  - first flux [mm/d]
%               In2  - second flux [mm/d]
%
% WK, 08/10/2018

out = max(In1-In2,0);

end

