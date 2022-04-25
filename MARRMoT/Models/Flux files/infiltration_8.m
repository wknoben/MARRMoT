function [out] = infiltration_8(S,Smax,fin)
%infiltration_8

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Infiltration into storage is equal the inflow when current 
%               storage is under the maximum storage, 
%               and zero when storage reaches maximum capacity 
% Constraints:  f <= fin
% @(Inputs):    
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]

out = (S < Smax) .* fin;

end