function [out] = smoothThreshold_temperature_logistic(T,Tt,r)
%smoothThreshold_temperature_logistic Logisitic smoother for temperature threshold functions.

% Copyright (C) 2018 Wouter J.M. Knoben
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%   Smooths the transition of threshold functions of the form:
%
%   Snowfall = { P, if T <  Tt
%              { 0, if T >= Tt
%
%   By transforming the equation above to Sf = f(P,T,Tt,r):
%   Sf = P * 1/ (1+exp((T-Tt)/r))
%
%   Inputs:
%   P       : current precipitation
%   T       : current temperature
%   Tt      : threshold temperature below which snowfall occurs
%   r       : [optional] smoothing parameter rho, default = 0.01
%
%   NOTE: this function only outputs the multiplier. This needs to be
%   applied to the proper flux utside of this function.

% Check for inputs and use defaults if not provided
% NOTE: this is not very elegant, but it is more than a factor 10 faster then: 
% if ~exist('r','var'); r = 0.01; end
if nargin == 2
    r = 0.01;
end

% Calculate multiplier
out = 1 ./ (1+exp((T-Tt)/(r)));

end

