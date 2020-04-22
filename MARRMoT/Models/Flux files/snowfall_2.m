function [func] = snowfall_2(~)
%snowfall_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Snowfall based on a temperature threshold interval
% Constraints:  -
% @(Inputs):    p1   - midpoint of the combined rain/snow interval [oC]
%               p2   - length of the mixed snow/rain interval [oC]
%               T    - current temperature [oC]
%               In   - incoming precipitation flux [mm/d]
%
% WK, 08/10/2018

func = @(In,T,p1,p2) min(In,max(0,In.*(p1+0.5*p2-T)/p2));

end

