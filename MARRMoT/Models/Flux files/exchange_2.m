function [out] = exchange_2(p1,S1,S1max,S2,S2max)
%exchange_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Water exchange based on relative storages
% Constraints:  -
% @(Inputs):    p1    - base exchange rate [mm/d]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]
%
% WK, 07/10/2018

out = p1*(S1/S1max - S2/S2max); 

end

