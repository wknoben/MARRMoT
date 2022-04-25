function [out] = saturation_9(In,S,St,varargin)
%saturation_9 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Deficit store: Saturation excess from a store that has 
%               reached maximum capacity
% Constraints:  -
% @(Inputs):    In   - incoming flux [mm/d]
%               S    - current storage [mm]
%               St   - threshold for flow generation [mm], 0 for deficit
%               store
%               varargin(1) - smoothing variable r (default 0.01)
%               varargin(2) - smoothing variable e (default 5.00)

if size(varargin,2) == 0
    out = In.*smoothThreshold_storage_logistic(S,St);
elseif size(varargin,2) == 1
    out = In.*smoothThreshold_storage_logistic(S,St,varargin(1));
elseif size(varargin,2) == 2
    out = In.*smoothThreshold_storage_logistic(S,St,varargin(1),varargin(2));    
end

end

