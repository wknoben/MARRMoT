function [out] = evap_11(S,Smax,Ep)
%evap_11 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation quadratically related to current soil moisture
% Constraints:  f >= 0
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%
% WK, 05/10/2018

out = max(0,(2*S/Smax-(S/Smax)^2)*Ep);

end

