function [func] = rainfall_1(varargin)
%rainfall_1 
% varargin(1): value of smoothing variable r (default 0.01)
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Rainfall based on temperature threshold
% Constraints:  -
% @(Inputs):    p1   - temperature threshold above which rainfall occurs [oC]
%               T    - current temperature [oC]
%               In   - incoming precipitation flux [mm/d]
%
% WK, 08/10/2018

if size(varargin,2) == 0
    func = @(In,T,p1) In.*(1-smoothThreshold_temperature_logistic(T,p1));
elseif size(varargin,2) == 1
    func = @(In,T,p1) In.*(1-smoothThreshold_temperature_logistic(T,p1,varargin(1)));    
end


end

