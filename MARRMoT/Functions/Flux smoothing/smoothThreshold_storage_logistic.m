function [out] = smoothThreshold_storage_logistic(S,Smax,r,e)
%smoothThreshold_storage_logistic Logisitic smoother for storage threshold functions.

% Copyright (C) 2018 Wouter J.M. Knoben
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%   Smooths the transition of threshold functions of the form:
%
%   Q = { P, if S = Smax
%       { 0, if S < Smax
%
%   By transforming the equation above to Q = f(P,S,Smax,e,r):
%   Q = P * 1/ (1+exp((S-Smax+r*e*Smax)/(r*Smax)))
%
%   Inputs:
%   S       : current storage
%   Smax    : maximum storage
%   r       : [optional] smoothing parameter rho, default = 0.01
%   e       : [optional] smoothing parameter e, default 5
%
%   NOTE: this function only outputs the multiplier. This needs to be
%   applied to the proper flux utside of this function.
%
%   NOTE: can be applied for temperature thresholds as well (i.e. snow
%   modules). This simply means that S becomes T, and Smax T0.

% Check for inputs and use defaults if not provided
% NOTE: this is not very elegant, but it is more than a factor 10 faster then: 
% if ~exist('r','var'); r = 0.01; end
% if ~exist('e','var'); e = 5.00; end
if nargin == 2
    r = 0.01;
    e = 5.00;
elseif nargin == 3
    e = 5.00;
end

% Calculate multiplier
Smax = max(Smax,0);   % this avoids numerical instabilities when Smax<0
if r*Smax == 0
    out = 1 ./ (1+exp((S-Smax+r*e*Smax)/(r)));
else
    out = 1 ./ (1+exp((S-Smax+r*e*Smax)/(r*Smax)));
end

end

