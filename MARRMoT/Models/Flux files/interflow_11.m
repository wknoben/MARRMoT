function [out] = interflow_11(p1,p2,S,dt,varargin)
%interflow_11 
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
% Description:  Constant interflow if storage exceeds a threshold
% Constraints:  f <= (S-p2_/dt
% @(Inputs):    p1   - interflow rate [mm/d]
%               p2   - storage threshold for flow generation [mm]
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 09/10/2018

if size(varargin,2) == 0
    out = min(p1,(S-p2)/dt).*(1-smoothThreshold_storage_logistic(S,p2));
elseif size(varargin,2) == 1
    out = min(p1,(S-p2)/dt).*(1-smoothThreshold_storage_logistic(S,p2,varargin(1)));
elseif size(varargin,2) == 2
    out = min(p1,(S-p2)/dt).*(1-smoothThreshold_storage_logistic(S,p2,varargin(1),varargin(2)));    
end

end

