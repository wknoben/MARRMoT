function [func] = baseflow_8(~)
%baseflow_8 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Exponential scaled outflow from a deficit store
% Constraints:  S <= Smax
% @(Inputs):    p1   - base outflow rate [mm/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%
% WK, 05/10/2018

func = @(p1,p2,S,Smax) p1.*(exp(p2.*min(1,max(S,0)./Smax))-1);

end

