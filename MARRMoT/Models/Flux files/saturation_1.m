function [out] = saturation_1(In,S,Smax,varargin)
%saturation_1 
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
% Description:  Saturation excess from a store that has reached maximum capacity
% Constraints:  -
% @(Inputs):    In   - incoming flux [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%
% WK, 09/10/2018

if size(varargin,2) == 0
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax));
elseif size(varargin,2) == 1
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1)));
elseif size(varargin,2) == 2
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1),varargin(2)));    
end

end

