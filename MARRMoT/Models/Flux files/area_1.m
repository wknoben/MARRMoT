function [out] = area_1(p1,p2,S,Smin,Smax,varargin)
%area_1 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Auxiliary function that calculates a variable contributing area.
% Constraints:  A <= 1
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smin - minimum contributing storage [mm]
%               Smax - maximum contributing storage [mm]
%               varargin(1) - smoothing variable r (default 0.01)
%               varargin(2) - smoothing variable e (default 5.00)

if size(varargin,2) == 0
    out = min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin));                       % default smoothing
elseif size(varargin,2) == 1
    out = min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin,varargin(1)));           % user-specified smoothing
elseif size(varargin,2) == 2
     out = min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
            (1-smoothThreshold_storage_logistic(S,Smin,varargin(1),varargin(2))); % user-specified smoothing
end

end

