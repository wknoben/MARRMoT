function [out] = exchange_2(p1,S1,S1max,S2,S2max)
%exchange_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Water exchange based on relative storages
% Constraints:  -
% @(Inputs):    p1    - base exchange rate [mm/d]
%               S1    - current storage in S1 [mm]
%               S1max - maximum storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]

out = p1*(S1/S1max - S2/S2max); 

end

