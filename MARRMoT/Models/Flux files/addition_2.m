function [func] = addition_2(~)
% addition_2
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Sums two fluxes
% Constraints:  None
% @(Inputs):    fin1   - incomming flux 1 [mm/d]
%               fin2   - incomming flux 2 [mm/d]

func = @(fin1,fin2) (fin1 + fin2);

end