function [out] = interflow_4(p1,p2,S)
%interflow_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Combined linear and scaled quadratic interflow
% Constraints:  f <= S
%               S >= 0     - prevents numerical issues with complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - scaling factor [mm-1 d-1]
%               S    - current storage [mm]

out = min(max(S,0),p1*max(S,0)+p2*max(S,0)^2);

end

