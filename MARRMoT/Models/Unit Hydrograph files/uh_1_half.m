function [ UH ] = uh_1_half( d_base, delta_t )
%uh_1_half Unit Hydrograph [days] with half a bell curve. GR4J-based

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
%   Unit hydrograph spreads the input volume over a time period x4.
%   Percentage of input returned only increases. 
%   I.e. d_base = 3.8 [days], delta_t = 1: 
%   UH(1) = 0.04  [% of inflow]
%   UH(2) = 0.17
%   UH(3) = 0.35
%   UH(4) = 0.45

%%TIME STEP SIZE
delay = d_base/delta_t;
if delay == 0; delay = 1; end       % any value below t = 1 means no delay,
                                    % but zero leads to problems
tt = 1:ceil(delay);                 % Time series for which we need UH 
                                    % ordinates [days]

%%EMPTIES
SH = zeros(1,length(tt)+1); SH(1) = 0;
UH = zeros(1,length(tt));

%%UNIT HYDROGRAPH
for t = tt
    if     t <  delay; SH(t+1) = (t./delay).^(5./2);
    elseif t >= delay; SH(t+1) = 1;
    end
    
    UH(t) = SH(t+1)-SH(t);
end

   UH(2,:) = zeros(size(UH));

end

