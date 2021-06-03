function [val,idx] = of_log_NSE(obs,sim,idx)
% of_log_NSE Calculates the Nash-Sutcliffe Efficiency (Nash & Sutcliffe, 1970) 
% of the log of simulated streamflow. Ignores time steps with negative flow
% values. Adds a constant e of 1/100 of mean(obs) to avoid issues with zero
% flows (Pushpalatha et al. 2012).
%
% Copyright (C) 2021 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
%
% Out:
% val       - objective function value          [1x1]
% idx       - optional vector of indices to use for calculation
%
% Nash, J. E.; Sutcliffe, J. V. (1970). "River flow forecasting through 
% conceptual models part I — A discussion of principles". Journal of 
% Hydrology. 10 (3): 282–290. Bibcode:1970JHyd...10..282N. 
% doi:10.1016/0022-1694(70)90255-6.
%
% Pushpalatha, R.; Perrin, C.; le Moine, N. and Andréassian V. (2012). "A
% review of efficiency criteria suitable for evaluating low-flow
% simulations". Journal of Hydrology. 420-421, 171-182. 
% doi:10.1016/j.jhydrol.2011.11.055

%% Check inputs and set defaults
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

%% filter to only selected indices
obs = obs(idx);
sim = sim(idx);                                            

%% Find the constant e
e = mean(obs)/100;

%% Apply constant and transform flows
obs = log(obs+e);
sim = log(sim+e);

%% Calculate metric
top = sum((sim - obs).^2);
bot = sum((obs - mean(obs)).^2);
val = 1 - (top/bot);
end

