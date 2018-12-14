function [func] = baseflow_9(~)
%baseflow_9 Linear flow above a threshold
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Linear outflow from a reservoir if a storage threshold is exceeded
% Constraints:  -
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - storage threshold for flow generation [mm]
%               S    - current storage [mm]
%
% WK, 05/10/2018

func = @(p1,p2,S) p1.*max(0,S-p2);

end

