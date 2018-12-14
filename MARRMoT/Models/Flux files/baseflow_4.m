function [func] = baseflow_4(~)
%baseflow_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Exponential outflow from deficit store
% Constraints:  - 
% @(Inputs):    p1   - base outflow rate [mm/d]
%               p2   - exponent parameter [mm-1]
%               S    - current storage [mm]
%
% WK, 05/10/2018

func = @(p1,p2,S) p1*exp(-1*p2*S);

end

