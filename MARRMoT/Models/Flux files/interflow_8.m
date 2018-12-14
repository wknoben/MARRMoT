function [func] = interflow_8(~)
%interflow_8 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Linear interflow if storage exceeds a threshold
% Constraints:  f = 0 for S < p2
% @(Inputs):    S    - current storage [mm]
%               p1   - time coefficient [d-1]
%               p2   - storage threshold before flow occurs [mm]
%
% WK, 09/10/2018

func = @(S,p1,p2)  max(0,p1*(S-p2));

end

