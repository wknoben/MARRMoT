function [func] = recharge_5(~)
%recharge_5 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Recharge to fulfil evaporation demand if the receiving 
%               store is below a threshold
% Constraints:  -
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - storage threshold in S2 [mm]
%               S1   - current storage in S1 [mm]
%               S2   - current storage in S1 [mm]
%
% WK, 08/10/2018

func = @(p1,p2,S1,S2) p1.*S1.*(1-min(1,S2./p2));

end

