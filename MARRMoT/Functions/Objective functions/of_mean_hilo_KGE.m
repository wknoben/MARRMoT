function [val,c,w] = of_mean_hilo_KGE(obs,sim,varargin)
% of_mean_hilo_KGE Calculates the average Kling-Gupta Efficiency of the 
% normal and inverse of simulated streamflow (Gupta et al, 2009, 
% Pushpalatha et al, 2012).
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% varargin  - optional weights of components    [3x1],[3x1] [high, low]
%
% Out:
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         [3x2] [high; low]
% w         - weights    [wr,wa,wb]             [3x2] [high; low]
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
% Decomposition of the mean squared error and NSE performance criteria: 
% Implications for improving hydrological modelling. Journal of Hydrology, 
% 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
%
% Pushpalatha, R., Perrin, C., Moine, N. Le, & Andréassian, V. (2012). A 
% review of efficiency criteria suitable for evaluating low-flow 
% simulations. Journal of Hydrology, 420–421, 171–182. 
% https://doi.org/10.1016/j.jhydrol.2011.11.055

%% call individual KGE functions (they have their own error checking)
if ~isempty(varargin)
    [val1,c1,w1] = of_KGE(obs,sim,varargin{1});
    [val2,c2,w2] = of_inverse_KGE(obs,sim,varargin{2});
else
    [val1,c1,w1] = of_KGE(obs,sim);
    [val2,c2,w2] = of_inverse_KGE(obs,sim);
end

%% calculate value
val = 0.5*(val1+val2);      % weighted KGE
c   = [c1; c2];             % components from high and low KGE
w   = [w1; w2];             % weights from high and low KGE

end