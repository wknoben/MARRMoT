function [func] = subtraction_2(~)
% subtraction_2
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Subtracts two fluxes
% Constraints:  None
% @(Inputs):    fin1   - incomming flux 1
%               fout1  - outgoing flux 1

func = @(fin1,fout1) (fin1 - fout1);

end