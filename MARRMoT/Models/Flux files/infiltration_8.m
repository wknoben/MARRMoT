function [func] = infiltration_8(~)
%infiltration_8
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Infiltration into storage is equal the inflow when current 
%               storage is under the maximum storage, 
%               and zero when storage reaches maximum capacity 
% Constraints:  f <= fin
% @(Inputs):    
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]

func = @(S,Smax,fin) ((S < Smax) .* fin);

end