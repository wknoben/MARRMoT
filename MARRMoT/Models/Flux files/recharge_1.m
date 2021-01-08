function [out] = recharge_1(p1,S,Smax,flux)
%recharge_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Recharge as scaled fraction of incoming flux
% Constraints:  -
% @(Inputs):    p1   - fraction of flux that is recharge [-]
%               S    - current storage [mm]
%               Smax - maximum contributing storage [mm]
%               flux - incoming flux [mm/d]
%
% WK, 08/10/2018

out = p1*S/Smax*flux;

end

