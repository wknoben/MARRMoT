function [out] = melt_1(p1,p2,T,S,dt)
%melt_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Snowmelt from degree-day-factor 
% Constraints:  f <= S/dt
% @(Inputs):    p1   - degree-day factor [mm/oC/d]
%               p2   - temperature threshold for snowmelt [oC]
%               T    - current temperature [oC]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 08/10/2018

out = max(min(p1*(T-p2),S/dt),0);

end

