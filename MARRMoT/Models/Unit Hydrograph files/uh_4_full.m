function [ UH ] = uh_4_full( d_base, delta_t )
%uh_4_half Unit Hydrograph [days] with a triangle (linear)

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
%   Percentage runoff goes up, peaks, and goes down again.
%   I.e. d_base = 3.8 [days], delta_t = 1: 
%   UH(1) = 0.14  [% of inflow]
%   UH(2) = 0.41
%   UH(3) = 0.36
%   UH(4) = 0.09

%%TIME STEP SIZE
delay = d_base/delta_t;
if delay == 0; delay = 1; end       % any value below t = 1 means no delay, 
                                    % but zero leads to problems
tt = 1:ceil(delay);                 % time series for which we need UH 
                                    % ordinates [days]

%%UNIT HYDROGRAPH
% The area under the unit hydrograph by definition sums to 1. Thus the area
% is S(t=0 to t = delay) t*[ff: fraction of flow to move per time step] dt
% Analytical solution is [1/2 * t^2 + c]*ff, with c = 0. 
% Here, we use two half triangles t make one big one, so the area of half a
% triangle is 0.5. Thus the fraction of flow step size is:
ff  = 0.5/(0.5*(0.5*delay)^2);           
d50 = 0.5*delay;

%%TRIANGLE FUNCTION
tri = @(t) max(ff.*(t-d50).*sign(d50-t)+ff.*d50,0);

%%EMPTIES
UH = zeros(1,length(tt));

%%UNIT HYDROGRAPH
for t = 1:length(tt) 
    UH(t) = integral(tri,t-1,t);
end

%%ENSURE UH SUMS TO 1
tmp_diff   = 1-sum(UH);
tmp_weight = UH./sum(UH);
UH         = UH + tmp_weight.*tmp_diff;

UH(2,:) = zeros(size(UH));

end

