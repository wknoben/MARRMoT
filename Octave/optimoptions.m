function [ opt ] = optimoptions( varargin )
%optimoptions Octave only! - Sets a simple setting for the root-finding
%method. An equivalent to Matlab's optimoptions is not (yet) available in
%Octave, and this function ensures the code still works.
%
% NOTE: this is NOT an Octave implementation of Matlab's optimoptions. DO
% NOT INCLUDE this file with MARRMoT's Matlab distribution because Matlab
% will preferentially use this file instead of its own built-in 
% optimoptions.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

opt = optimset('MaxFunEvals',1000);  


end

