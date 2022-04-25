function [out] = soilmoisture_1(S1,S1max,S2,S2max)
%soilmoisture_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Water rebalance to equal relative storage (2 stores)
% Constraints:  -
% @(Inputs):    S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]

out = ((S2.*S1max-S1.*S2max)/(S1max+S2max)).*smoothThreshold_storage_logistic(S1./S1max,S2./S2max);

end

