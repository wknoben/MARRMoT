function [out] = interception_3(p1,In)
%interception_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Interception excess after a fraction is intercepted
% Constraints:  -
% @(Inputs):    p1   - fraction throughfall [-]
%               In   - incoming flux [mm/d]

out = p1.*In;

end

