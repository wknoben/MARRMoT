function [func] = evap_18(~)
%evap_18 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Exponentially declining evaporation from deficit store
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - linear scaling parameter [-]
%               p3   - storage scaling parameter [mm]
%               S    - current storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%
% WK, 05/10/2018

func = @(p1,p2,p3,S,Ep) p1.*exp(-1.*p2.*S./p3).*Ep;

end

