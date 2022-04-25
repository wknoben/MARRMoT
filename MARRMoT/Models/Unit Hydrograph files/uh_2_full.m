function [ UH ] = uh_2_full( d_base,delta_t )
%uh_2_full Unit Hydrograph [days] with a full bell curve. GR4J-based

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
%   Unit hydrograph spreads the input volume over a time period 2*x4.
%   Percentage of input returned goes up (till x4), then down again.
%   I.e. d_base = 3.8 [days], delta_t = 1:
%   UH(1) = 0.02  [% of inflow]
%   UH(2) = 0.08
%   UH(3) = 0.18
%   UH(4) = 0.29
%   UH(5) = 0.24
%   UH(6) = 0.14
%   UH(7) = 0.05
%   UH(8) = 0.00

%%TIME STEP SIZE
delay = d_base/delta_t;
tt = 1:2*ceil(delay); % time series for which we need UH ordinates [days]

%%EMPTIES
SH = zeros(1,length(tt)+1); SH(1) = 0;
UH = zeros(1,length(tt));

%%UNIT HYDROGRAPH
for t = tt
    if     (t <= delay)
        SH(t+1) = 0.5*(t./delay).^(5./2);
    elseif (t > delay) && (t < 2*delay); 
        SH(t+1) = 1-0.5*(2-t./delay).^(5./2);
    elseif (t >= 2*delay);
        SH(t+1) = 1;
    end
    
    UH(t) = SH(t+1)-SH(t);
end

    UH(2,:) = zeros(size(UH));

end

