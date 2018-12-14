function [func] = infiltration_5(~)
%infiltration_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Maximum infiltration rate non-linearly based on relative deficit and storage
% Constraints:  S2 >= 0     - prevents complex numbers
%               f <= 10^9   - prevents numerical issues with Inf outcomes
% @(Inputs):    p1    - base infiltration rate [mm/d]
%               p2    - exponential scaling parameter [-]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]
%
% WK, 07/10/2018

func = @(p1,p2,S1,S1max,S2,S2max) max(0,min(10^9,p1.*(1-S1./S1max).*max(0,S2./S2max).^(-1.*p2)));

end

