function [func] = infiltration_2(~)
%infiltration_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  0 <= f <= S2/dt
% @(Inputs):    p1    - maximum infiltration rate [mm,/d]
%               p2    - exponential scaling parameter [-]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               flux  - reduction of infiltration rate by infiltration 
%                       demand already fulfilled elsewhere [mm/d]
%               S2    - storage available for infiltration [mm]
%               dt    - time step size [d]
%
% WK, 07/10/2018

func = @(p1,p2,S1,S1max,flux,S2,dt) max(min((p1.*exp(-1*p2*S1./S1max))-flux,S2/dt),0);

end

