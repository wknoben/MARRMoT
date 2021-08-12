function [out] = interception_2(In,p1)
%interception_2 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Interception excess after a constant amount is intercepted
% Constraints:  f >= 0
% @(Inputs):    In   - incoming flux [mm/d]
%               p1   - interception and evaporation capacity [mm/d]

out = max(In-p1,0);

end

