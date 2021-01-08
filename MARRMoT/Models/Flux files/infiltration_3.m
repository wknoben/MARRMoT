function [out] = infiltration_3(In,S,Smax,varargin)
%infiltration_3 
% varargin(1): value of smoothing variable r (default 0.01)
% varargin(2): value of smoothing variable e (default 5.00)
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Infiltration to soil moisture of liquid water stored in snow pack
% Constraints:  -
% @(Inputs):    In   - incoming flux [mm/d]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%
% WK, 07/10/2018

if size(varargin,2) == 0
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax));
elseif size(varargin,2) == 1
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1)));
elseif size(varargin,2) == 2
    out = In.*(1-smoothThreshold_storage_logistic(S,Smax,varargin(1),varargin(2)));    
end

end

