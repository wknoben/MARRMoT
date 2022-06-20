function [val,c,idx,w] = of_mean_hilo_root5_KGE(obs,sim,idx,w)
% of_mean_hilo_KGE Calculates the average Kling-Gupta Efficiency (Gupta et
% al., 2009) of the untransformed and fifth root (Chiew et al., 1993) of
% streamflow. Ingores timesteps with negative flow values.

% Copyright (C) 2021 L. Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
% w         - optional weights of components    {[3x1] [3x1]} {high low}
                                            %or [3x1] high = low
%
% Out:
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         {[3x1] [3x1]} {high low}
% idx       - indices used for the calculation
% w         - weights    [wr,wa,wb]             {[3x1] [3x1]} {high low}
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
% Decomposition of the mean squared error and NSE performance criteria: 
% Implications for improving hydrological modelling. Journal of Hydrology, 
% 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
%
% Chiew, F.H.S., Stewardson, M.J., McMahon, T.A., 1993. Comparison of six
% rainfall-runoff modelling approaches. J. Hydrol.
% https://doi.org/10.1016/0022-1694(93)90073-I

%% Check inputs and select timesteps
if nargin < 2
    error('Not enugh input arguments')    
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);

%% Set weights
if nargin < 4 || isempty(w)
    w1 = [1 1 1];
    w2 = [1 1 1];
else
    if iscell(w) && size(w) == 2
        w1 = w{1}; w2 = w{2};
    elseif isarray(w)
        w1 = w; w2 = w;
    else
        error(['w should be either'...
               'a cell array of size 2, or'...
               'an array of size [1,3]']);
    end
end

%% call individual KGE functions (they have their own error checking)
[val1,c1,~,w1] = of_KGE(obs,sim,[],w1);
[val2,c2,~,w2] = of_root5_KGE(obs,sim,[],w2);

%% calculate value
val = 0.5*(val1+val2);      % weighted KGE
c   = {c1 c2};              % components from high and low KGE
w   = {w1 w2};              % weights from high and low KGE

end