function [func] = excess_1(~)
%excess_1
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Storage excess when store size changes (returns flux [mm/d])
% Constraints:  f >= 0
% @(Inputs):    So   - 'old' storage [mm]
%               Smax - 'new' maximum storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(So,Smax,dt) max((So-Smax)/dt,0);

end

