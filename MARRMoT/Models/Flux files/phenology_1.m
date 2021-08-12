function [ out ] = phenology_1(T,p1,p2,Ep)
%phenology_1 Corrects Ep for phenology effects.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Phenology-based correction factor for potential
%               evapotranspiration (returns flux [mm/d])
% Constraints:  -
% @(Inputs):    T    - current temperature [oC]
%               p1   - temperature threshold where evaporation stops [oC]
%               p2   - temperature threshold above which corrected Ep = Ep
%               Ep   - current potential evapotranspiration [mm/d]

out = min(1,max(0,(T-p1)/(p2-p1)))*Ep;

end

