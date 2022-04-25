function [ UH ] = uh_8_delay( t_delay, delta_t )
%uh_8_delay Unit Hydrograph [days] with a pure delay (no transformation).

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%   Inputs
%   t_delay - flow delay [d]
%   delta_t - time step size [d] 
%
%   Output
%   UH      - unit hydrograph [nx2]
%               uh's first row contains coeficients to splut flow at each
%               of n timesteps forward, the second row contains zeros now,
%               these are the still-to-flow values.
%
%   Unit hydrograph shifts the input volume over a time period.
%   Input is spread over maximum 2 time steps.
%   I.e. t_delay = 3.8 [days], delta_t = 1:
%   UH(1) = 0.00  [% of inflow]
%   UH(2) = 0.00
%   UH(3) = 0.00
%   UH(4) = 0.20
%   UH(5) = 0.80

%%TIME STEP SIZE
delay = t_delay/delta_t;

%%UNIT HYDROGRAPH
% The input is only shifted in time, not transformed, so we only need two
% ordinates:
ord1 = 1-t_delay+floor(t_delay);
ord2 = t_delay-floor(t_delay);

% Flow appears from this time step (t+t_start; a delay of 1 time step means
% flow doesn't appear on t=1, but starts on t=1+1):
t_start = floor(delay);

% Unit Hydrograph
UH(1+t_start) = ord1;
UH(1+t_start+1) = ord2;

UH(2,:) = zeros(size(UH));

end

