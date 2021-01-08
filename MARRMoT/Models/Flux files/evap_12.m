function [out] = evap_12(S,p1,Ep)
%evap_11 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation from deficit store, with exponential decline as
%               deficit goes below a threshold
% Constraints:  -
% @(Inputs):    S    - current storage [mm]
%               p1   - wilting point [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%
% WK, 05/10/2018

out = min(1,exp(2*(1-S/p1)))*Ep;

end

