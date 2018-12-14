function [func] = area_1(varargin)
%area_1 
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
% Description:  Auxiliary function that calculates a variable contributing area.
% Constraints:  A <= 1
% @(Inputs):    p1   - linear scaling parameter [-]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smin - minimum contributing storage [mm]
%               Smax - maximum contributing storage [mm]
%
% WK, 09/10/2018

% Return a function based on the provided varargin
if size(varargin,2) == 0
    func = @(p1,p2,S,Smin,Smax) ...
        min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin));                       % default smoothing
elseif size(varargin,2) == 1
    func = @(p1,p2,S,Smin,Smax) ...
        min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin,varargin(1)));           % user-specified smoothing
elseif size(varargin,2) == 2
     func = @(p1,p2,S,Smin,Smax) ...
        min(1,p1.*(max(0,S-Smin)./(Smax-Smin)).^p2).*...
        (1-smoothThreshold_storage_logistic(S,Smin,varargin(1),varargin(2))); % user-specified smoothing
end

end

