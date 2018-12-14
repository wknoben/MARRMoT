function [func] = saturation_11(varargin)
%saturation_11 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% varargin(1): value of smoothing variable r (default 0.01)
% varargin(2): value of smoothing variable e (default 5.00)
%
% Anonymous function
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
%
% WK, 09/10/2018

if size(varargin,2) == 0
    func = @(p1,p2,S,Smin,Smax,In) ...
        In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin));
elseif size(varargin,2) == 1
    func = @(p1,p2,S,Smin,Smax,In) ...
        In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin,varargin(1)));
elseif size(varargin,2) == 2
     func = @(p1,p2,S,Smin,Smax,In) ...
        In.*min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin,varargin(1),varargin(2)));   
end


end

