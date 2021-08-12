function [out] = saturation_11(p1,p2,S,Smin,Smax,In,varargin)
%saturation_11 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (min exponential variant)
% Constraints:  f <= In
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smin - minimum contributing storage [mm]
%               Smax - maximum contributing storage [mm]
%               In   - incoming flux [mm/d]
%               varargin(1) - smoothing variable r (default 0.01)
%               varargin(2) - smoothing variable e (default 5.00)

if size(varargin,2) == 0
    out = In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin));
elseif size(varargin,2) == 1
    out = In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin,varargin(1)));
elseif size(varargin,2) == 2
     out = In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin,varargin(1),varargin(2)));   
end


end

