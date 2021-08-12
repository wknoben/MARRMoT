function [ UH ] = uh_7_uniform( d_base, delta_t )
%uh_7_uniform Unit Hydrograph [days] with uniform spread

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%   Inputs
%   d_base  - time base of routing delay [d]
%   delta_t - time step size [d]    
%
%   Output
%   UH      - unit hydrograph [nx2]
%               uh's first row contains coeficients to splut flow at each
%               of n timesteps forward, the second row contains zeros now,
%               these are the still-to-flow values.
%
%   Unit hydrograph spreads the input volume over a time period delay.
%   I.e. d_base = 3.8 [days], delta_t = 1:  
%   UH(1) = 0.26  [% of inflow]
%   UH(2) = 0.26
%   UH(3) = 0.26
%   UH(4) = 0.22

%%TIME STEP SIZE
delay = d_base/delta_t;
tt = 1:ceil(delay); % time series for which we need UH ordinates [days]

%%EMPTIES
UH = NaN.*zeros(1,length(tt));

%%FRACTION FLOW
ff = 1/delay; % fraction of flow per time step

%%UNIT HYDROGRAPH
for t=1:ceil(delay)
    if t < delay
        UH(t) = ff;
    else
        UH(t) = mod(delay,t-1)*ff;
    end
end

UH(2,:) = zeros(size(UH));

end

