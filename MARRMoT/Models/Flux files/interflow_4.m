function [out] = interflow_4(p1,p2,S)
%interflow_4 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Combined linear and scaled quadratic interflow
% Constraints:  f <= S
%               S >= 0     - prevents numerical issues with complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - scaling factor [mm-1 d-1]
%               S    - current storage [mm]
%
% WK, 08/10/2018

out = min(max(S,0),p1*max(S,0)+p2*max(S,0)^2);

end

