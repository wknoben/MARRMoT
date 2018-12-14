function [func] = saturation_9(varargin)
%saturation_9 
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
% Description:  Deficit store: Saturation excess from a store that has 
%               reached maximum capacity
% Constraints:  -
% @(Inputs):    In   - incoming flux [mm/d]
%               S    - current storage [mm]
%               St   - threshold for flow generation [mm], 0 for deficit
%               store
%
% WK, 09/10/2018

if size(varargin,2) == 0
    func = @(In,S,St) In.*smoothThreshold_storage_logistic(S,St);
elseif size(varargin,2) == 1
    func = @(In,S,St) In.*smoothThreshold_storage_logistic(S,St,varargin(1));
elseif size(varargin,2) == 2
    func = @(In,S,St) In.*smoothThreshold_storage_logistic(S,St,varargin(1),varargin(2));    
end

end

