function [out] = infiltration_7(p1,p2,S,Smax,In,varargin)
%infiltration_7: infiltration_1 with customization

% Copyright (C) 2021 Clara Brandes, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  f <= fin
% @(Inputs):    p1   - maximum infiltration rate [mm,/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]

pre_smoother = (min(p1.*exp((-1*p2*S)./Smax),In));
if size(varargin,2) == 0
    out = pre_smoother.*(1-smoothThreshold_storage_logistic(S,Smax));
elseif size(varargin,2) == 1
    out = pre_smoother.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1)));
elseif size(varargin,2) == 2
    out = pre_smoother.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1),varargin(2)));    
end

end