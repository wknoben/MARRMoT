function [func] = abstraction_1(~)
%abstraction_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Function for constant groundwater abstraction
% Constraints:  None, taken from a store with possible negative depth
% @(Inputs):    p1 - Abstraction rate [mm/d]
%
% WK, 05/10/2018

func = @(p1) p1;

end

