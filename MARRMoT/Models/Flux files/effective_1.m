function [out] = effective_1(In1,In2)
%effective_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  General effective flow (returns flux [mm/d])
% Constraints:  In1 > In2
% @(Inputs):    In1  - first flux [mm/d]
%               In2  - second flux [mm/d]

out = max(In1-In2,0);

end

