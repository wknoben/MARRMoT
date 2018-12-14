function [ out,UH,frac_routing_beyond_time_series ] = ...
                                        uh_6_gamma( in,n,k,t_end,delta_t )
%uh_6_gamma Unit Hydrograph [days] from gamma function.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
%   Inputs
%   n       = shape parameter [-]
%   k       = time delay for flow reduction by a factor e [d]
%   t_end   = length of time series [d]
%   delta_t = time step size [d]
%
%   Unit hydrograph spreads the input volume over a time period delay.
%   Percentage of input returned only decreases. 
%   I.e. n = 1, k = 3.8 [days], delta_t = 1:
%   UH(1) = 0.928  [% of inflow]
%   UH(2) = 0.067
%   UH(3) = 0.005
%   UH(4) = 0.000

%%INPUTS
if any(size(in)) > 1; error('UH input should be a single value.'); end

%%TIME STEP SIZE
tmax = t_end/delta_t;
tt   = 1:tmax;          % time series for which we need UH ordinates [days] 

%%EMPTIES
UH_full  = zeros(1,length(tt));
frac_routing_beyond_time_series = 0;

%%UNIT HYDROGRAPH
% The Unit Hydrograph follows a gamma distribution. For a given 
% delay time, the fraction of flow per time step is thus the integral of 
% t-1 to t of the gamma distrubtion. The curve has range [0,Inf>. 
% We need to choose a point at which to cap the integration, but this
% depends on the parameters n & k, and the total time step. We choose the 
% cutoff point at the time step where less than 0.1% of the peak flow 
% is still on route.

%%Unit hydrograph
for t = 1:length(tt)
    UH_full(t) = integral(@(x) 1./(k.*gamma(n)).*(x./k).^(n-1).* ...
                                    exp(-1.*x./k),(t-1)*delta_t,t*delta_t);
end

%%Find cutoff point where less than 0.1% of the peak flow is being routed
[max_val,max_here] = max(UH_full);
end_here = find(UH_full(max_here:end)./max_val<0.001,1) + max_here;

%%Take action depending on whether the distribution function exceeds the
%%time limit or not
if ~isempty(end_here) 
    %%Construct the Unit Hydrograph
    UH = UH_full(1:end_here);

    %%Account for the truncated part of the full UH.
    % find probability mass to the right of the cut-off point
    tmp_excess = 1-sum(UH);                                         

    % find relative size of each time step
    tmp_weight = UH_full(1:end_here)./sum(UH_full(1:end_here));     

    % distribute truncated probability mass proportionally to all elements 
    % of the routing vector
    UH = UH+tmp_weight.*tmp_excess;                                 
    
else
    %%Construct the Unit Hydrograph
    UH = UH_full;

    %%The UH is longer than the provided time series length. Track the
    %%percentage of flow that is routed beyond the simulation duration
    frac_routing_beyond_time_series = 1-sum(UH);
    
end

%%DISPERSE VOLUME
out = in.*UH;

end

