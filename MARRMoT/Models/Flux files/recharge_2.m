function [func] = recharge_2(~)
%recharge_2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Recharge as non-linear scaling of incoming flux
% Constraints:  S >= 0
% @(Inputs):    p1   - recharge scaling non-linearity [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               flux - incoming flux [mm/d]
%
% WK, 08/10/2018

func = @(p1,S,Smax,flux) flux*((max(S,0)/Smax)^p1);

end

