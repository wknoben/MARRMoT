function [func] = soilmoisture_1(~)
%soilmoisture_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Water rebalance to equal relative storage (2 stores)
% Constraints:  -
% @(Inputs):    S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]
%
% WK, 09/10/2018

func = @(S1,S1max,S2,S2max) ((S2.*S1max-S1.*S2max)/(S1max+S2max)).*smoothThreshold_storage_logistic(S1./S1max,S2./S2max);

end

