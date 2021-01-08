function [out] = saturation_8(p1,p2,S,Smax,In)
%saturation_8 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Saturation excess flow from a store with different degrees 
%               of saturation (min-max linear variant)
% Constraints:  -
% @(Inputs):    p1   - minimum fraction contributing area [-]
%               p2   - maximum fraction contributing area [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               In   - incoming flux [mm/d]
%
% WK, 09/10/2018

out = (p1+(p2-p1)*S/Smax)*In;

end

