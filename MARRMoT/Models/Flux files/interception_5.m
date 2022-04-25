function [out] = interception_5(p1,p2,In)
%interception_5 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Interception excess after a combined absolute amount and fraction are intercepted
% Constraints:  f >= 0
% @(Inputs):    p1   - fraction that is not throughfall [-]
%               p2   - constnat interception and evaporation [mm/d]
%               In   - incoming flux [mm/d]

out = max(p1.*In-p2,0);


end

