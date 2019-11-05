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

% Set options
opt = optimset('MaxFunEvals',1000);  

% Display a warning that Octave is used. User should disable this if they are intentionally using Octave
%disp('---')
%disp('You are currently using MARRMoT''s Octave version.') 
%disp('If this is intentional, please navigate to ./MARRMoT/Octave/optimoptions.m and %disable this warning on lines 19-23 of the file.')
%disp('If you want to use MARRMoT''s Matlab version, you must REMOVE the folder ./MARRMoT/Octave from your system. If you do not remove this folder, Matlab will not utilize its more powerful solver settings because it will preferentially navigate to the Octave folder and use those settings instead.')
%disp('---')
%pause

end

