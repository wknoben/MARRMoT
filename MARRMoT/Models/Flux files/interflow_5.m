function [out] = interflow_5(p1,S)
%interflow_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Linear interflow
% Constraints:  -
% @(Inputs):    p1   - time coefficient [d-1]
%               S    - current storage [mm]
%
% WK, 08/10/2018

out = p1.*S;

end

