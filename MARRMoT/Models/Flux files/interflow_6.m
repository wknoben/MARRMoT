function [out] = interflow_6(p1,p2,S1,S2,S2max)
%interflow_6 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Scaled linear interflow if a storage in the receiving store exceeds a threshold
% Constraints:  S2/S2max <= 1
% @(Inputs):    p1    - time coefficient [d-1]
%               p2    - threshold as fractional storage in S2 [-]
%               S1    - current storage in S1 [mm]
%               S2    - current storage in S2 [mm]
%               S2max - maximum storage in S2 [mm]

out = p1.*S1.*(min(1,S2./S2max)-p2)./(1-p2).*...
                            (1-smoothThreshold_storage_logistic(S2./S2max,p2));

end

