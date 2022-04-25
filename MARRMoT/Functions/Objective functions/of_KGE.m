function [val,c,idx,w] = of_KGE(obs, sim, idx, w)
% of_KGE Calculates Kling-Gupta Efficiency of simulated streamflow (Gupta
% et al, 2009). Ignores time steps with negative flow values.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
% w         - optional weights of components    [3x1]
%
% Out:
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         [3x1]
% idx       - indices used for the calculation
% w         - weights    [wr,wa,wb]             [3x1]
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
% Decomposition of the mean squared error and NSE performance criteria: 
% Implications for improving hydrological modelling. Journal of Hydrology, 
% 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003

%% Check inputs and select timesteps
if nargin < 2
    error('Not enugh input arguments')    
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);

%% Set weights
w_default = [1,1,1];          % default weights

% update defaults weights if needed  
if nargin < 4 || isempty(w)
    w = w_default;
else
    if ~min(size(w)) == 1 || ~max(size(w)) == 3                            % check weights variable for size
        error('Weights should be a 3x1 or 1x3 vector.')                    % or throw error        
    end
end   

%% calculate components
c(1) = corr(obs,sim);                                             % r: linear correlation
c(2) = std(sim)/std(obs);                                         % alpha: ratio of standard deviations
c(3) = mean(sim)/mean(obs);                                       % beta: bias 

%% calculate value
val = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2);    % weighted KGE

end