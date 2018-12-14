function [func] = baseflow_3(~)
%baseflow_3 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Empirical non-linear outflow from a reservoir
% Constraints:  None specified
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%
% WK, 05/10/2018

func = @(S,Smax) Smax^(-4)/4*(S^5);

end

