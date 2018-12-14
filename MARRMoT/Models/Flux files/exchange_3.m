function [func] = exchange_3(~)
%exchange_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Water exchange with infinite size store based on threshold
% Constraints:  -
% @(Inputs):    p1   - base leakage time delay [d-1]
%               p2   - threshold for flow reversal [mm]
%               S    - current storage [mm]
%
% WK, 07/10/2018

func = @(p1,S,p2) p1*(S-p2);

end

