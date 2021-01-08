function [out] = evap_20(p1,p2,S,Smax,Ep,dt)
%evap_20 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation limited by a maximum evaporation rate and 
%               scaled below a wilting point
% Constraints:  f <= Ep
%               f <= S/dt
% @(Inputs):    p1   - maximum evaporation rate [mm/d]
%               p2   - wilting point as fraction of Smax [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]
%
% WK, 05/10/2018

out = min([p1.*S./(p2.*Smax),Ep,S/dt]);

end

