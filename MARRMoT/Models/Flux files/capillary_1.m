function [func] = capillary_1(~)
%capillary_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Capillary rise: based on deficit in higher reservoir
% Constraints:  f <= S2/dt
% @(Inputs):    p1   - maximum capillary rise rate  [mm/d]
%               S1   - current storage in receiving store [mm]
%               S1max- maximum storage in receiving store [mm]
%               S2   - current storage in providing store [mm]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,S1,S1max,S2,dt) min(p1.*(1-S1/S1max),S2/dt);

end

