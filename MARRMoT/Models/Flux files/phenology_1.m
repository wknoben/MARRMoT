function [ func ] = phenology_1( ~ )
%phenology_1 Corrects Ep for phenology effects.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Phenology-based correction factor for potential
%               evapotranspiration (returns flux [mm/d])
% Constraints:  -
% @(Inputs):    T    - current temperature [oC]
%               p1   - temperature threshold where evaporation stops [oC]
%               p2   - temperature threshold above which corrected Ep = Ep
%               Ep   - current potential evapotranspiration [mm/d]
%
% WK, 08/10/2018

func = @(T,p1,p2,Ep) min(1,max(0,(T-p1)/(p2-p1)))*Ep;

end

