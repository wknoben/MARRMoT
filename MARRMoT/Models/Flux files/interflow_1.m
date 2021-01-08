function [out] = interflow_1(p1,S,Smax,flux)
%interflow_1 Creates function for interflow: based on incoming flux size
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Interflow as a scaled fraction of an incoming flux
% Constraints:  -
% @(Inputs):    p1   - linear scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               flux - incoming flux [mm/d]
%
% WK, 08/10/2018

out = p1*S/Smax*flux;

end

