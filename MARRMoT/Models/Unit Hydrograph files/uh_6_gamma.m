function [ UH ] = uh_6_gamma(n,k,delta_t)
%uh_6_gamma Unit Hydrograph [days] from gamma function.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%   Inputs
%   n       = shape parameter [-]
%   k       = time delay for flow reduction by a factor e [d]
%   delta_t = time step size [d]
%
%   Output
%   UH      - unit hydrograph [nx2]
%               uh's first row contains coeficients to splut flow at each
%               of n timesteps forward, the second row contains zeros now,
%               these are the still-to-flow values.
%
%   Unit hydrograph spreads the input volume over a time period delay.
%   Percentage of input returned only decreases. 
%   I.e. n = 1, k = 3.8 [days], delta_t = 1:
%   UH(1) = 0.928  [% of inflow]
%   UH(2) = 0.067
%   UH(3) = 0.005
%   UH(4) = 0.000

%%TIME STEP SIZE
%tmax = t_end/delta_t;
%tt   = 1:tmax;          % time series for which we need UH ordinates [days] 

%%EMPTIES
%UH_full2  = zeros(1,length(tt));
%frac_routing_beyond_time_series = 0;

%%UNIT HYDROGRAPH
% The Unit Hydrograph follows a gamma distribution. For a given 
% delay time, the fraction of flow per time step is thus the integral of 
% t-1 to t of the gamma distrubtion. The curve has range [0,Inf>. 
% We need to choose a point at which to cap the integration, but this
% depends on the parameters n & k, and the total time step. We choose the 
% cutoff point at the time step where less than 0.1% of the peak flow 
% is still on route.

t = 1;
while true
    % calculate the pdf of the gamma distr at this timestep
    UH(t) = integral(@(x) 1./(k.*gamma(n)).*(x./k).^(n-1).* ...
                                    exp(-1.*x./k),(t-1)*delta_t,t*delta_t);
    
    % if the new value of the UH is less than 0.1% of the peak, end the
    % loop. NOTE: this works because the gamma distr is monomodal, hence on
    % the way to the peak UH(t) = max(UH) > max(UH) * .001.
    if UH(t) < (max(UH) * .001); break; end
    
    % go to the next step
    t = t+1;
end

%%Account for the truncated part of the UH.
% find probability mass to the right of the cut-off point
tmp_excess = 1-sum(UH);                                         

% find relative size of each time step
tmp_weight = UH./sum(UH);     

% distribute truncated probability mass proportionally to all elements 
% of the routing vector
UH = UH+tmp_weight.*tmp_excess;

UH(2,:) = zeros(size(UH));
end

