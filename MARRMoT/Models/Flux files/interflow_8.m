function [out] = interflow_8(S,p1,p2)
%interflow_8 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Linear interflow if storage exceeds a threshold
% Constraints:  f = 0 for S < p2
% @(Inputs):    S    - current storage [mm]
%               p1   - time coefficient [d-1]
%               p2   - storage threshold before flow occurs [mm]

out =  max(0,p1*(S-p2));

end

