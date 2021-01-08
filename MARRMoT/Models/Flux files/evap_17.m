function [out] = evap_17(p1,S,Ep)
%evap_17 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Scaled evaporation from a store that allows negative values
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [mm-1]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%
% WK, 05/10/2018

out = 1./(1+exp(-1.*p1.*S)).*Ep;

end

