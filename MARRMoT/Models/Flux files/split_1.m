function [out] = split_1(p1,In)
%split_1 flow splitting

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Split flow (returns flux [mm/d])
% Constraints:  -
% @(Inputs):    p1   - fraction of flux to be diverted [-]
%               In   - incoming flux [mm/d]

out = p1.*In;

end

