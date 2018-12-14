function [func] = interflow_10(~)
%interflow_10 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Scaled linear interflow if storage exceeds a threshold
% Constraints:  f = 0, for S < p2
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - threshold for flow generation [mm]
%               p3   - linear scaling parameter [-]
%               S    - current storage [mm]
%
% WK, 09/10/2018

func = @(S,p1,p2,p3)  p1*max(0,S-p2)/(p3);

end

