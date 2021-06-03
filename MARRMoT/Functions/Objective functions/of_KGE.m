function [val,c,idx,w] = of_KGE(obs, sim, idx, w)
% of_KGE Calculates Kling-Gupta Efficiency of simulated streamflow (Gupta
% et al, 2009). Ignores time steps with -999 values.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
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

%% check inputs and set defaults
if nargin < 2
    error('Not enugh input arguments')
elseif nargin > 4
    error('Too many inputs.')    
end

% make sure inputs are vertical and have the same size
obs = obs(:);
sim = sim(:);
if ~size(obs) == size(sim)
    error('Time series not of equal size.')
end

% defaults
w_default = [1,1,1];          % weights
idx_exists = find(obs >= 0);  % time steps to use in calculating of value
% -999 is opten used to denote missing values in observed data. Therefore
% we check for all negative values, and ignore those. 

% update default indices if needed
if nargin < 3 || isempty(idx)
    idx = idx_exists;
else 
    idx = idx(:);
    if islogical(idx) && all(size(idx) == size(obs))
        idx = intersect(find(idx), idx_exists);
    elseif isnumeric(idx)
        idx = intersect(idx, idx_exists);
    else
        error(['Indices should be either ' ...
                'a logical vector of the same size of Qsim and Qobs, or '...
                'a numeric vector of indices']);
    end                                                      % use all non missing Q if idx is not provided otherwise
end


% update defaults weights if needed  
if nargin < 4 || isempty(w)
    w = w_default;
else
    if ~min(size(w)) == 1 || ~max(size(w)) == 3                            % check weights variable for size
        error('Weights should be a 3x1 or 1x3 vector.')                    % or throw error        
    end                                                         % use dafult weight is w is not provided
end   
%% filter to only selected indices
obs = obs(idx);
sim = sim(idx);                                            

%% calculate components
c(1) = corr(obs(idx),sim(idx));                                             % r: linear correlation
c(2) = std(sim(idx))/std(obs(idx));                                         % alpha: ratio of standard deviations
c(3) = mean(sim(idx))/mean(obs(idx));                                       % beta: bias 

%% calculate value
val = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2);    % weighted KGE

end