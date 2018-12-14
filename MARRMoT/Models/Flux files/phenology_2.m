function [func] = phenology_2(~)
%phenology_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Phenology-based maximum interception capacity (returns store size [mm])
% Constraints:  Implicit assumption: 0 <= p2 <= 1
% @(Inputs):    p1   - mean interception capacity [mm]
%               p2   - seasonal change as fraction of the mean [-]
%               p3   - time of maximum store size [d]
%               t    - current time step [-]
%               tmax - seasonal length [d]
%               dt   - time step size [d]
%
% WK, 08/10/2018

func = @(p1,p2,p3,t,tmax,dt) p1*(1+p2*sin(2*pi*(t*dt-p3)/tmax));

end

