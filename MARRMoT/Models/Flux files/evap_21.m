function [func] = evap_21(~)
%evap_21 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Threshold-based evaporation with constant minimum rate
% Constraints:  f <= S/dt
% @(Inputs):    p1   - wilting point (1st threshold) [mm]
%               p2   - 2nd threshold as fraction of wilting point [-]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

func = @(p1,p2,S,Ep,dt) min(max(p2,min(S./p1,1)).*Ep,S);

end

