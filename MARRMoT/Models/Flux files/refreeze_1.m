function [out] = refreeze_1(p1,p2,p3,T,S,dt)
%refreeze_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Refreezing of stored melted snow
% Constraints:  f <= S/dt
% @(Inputs):    p1   - reduction fraction of degree-day-factor [-]
%               p2   - degree-day-factor [mm/oC/d]
%               p3   - temperature threshold for snowmelt [oC]
%               T    - current temperature [oC]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = max(min(p1*p2*(p3-T), S/dt), 0);

end

