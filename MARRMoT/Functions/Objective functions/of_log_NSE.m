function [val] = of_log_NSE(obs,sim,varargin)
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
% varargin  - number of timesteps for warmup    [1x1]
%
% Out:
% val       - objective function value          [1x1]
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
elseif nargin > 3
    error('Too many inputs.')    
end

% Defaults
warmup = 0; % time steps to ignore when calculating

% update default warmup period if needed
if nargin == 3
    if size(varargin{1}) == [1,1]
        warmup = varargin{1};
    else
        error('Warm up period should be 1x1 scalar.')
    end
end

% check time series size and rotate one if needed
if checkTimeseriesSize(obs,sim) == 0
    error('Time series not of equal size.')
    
elseif checkTimeseriesSize(obs,sim) == 2
    sim = sim';                                                             % 2 indicates that obs and sim are the same size but have different orientations
end

%% Apply warmup period
obs = obs(1+warmup:end);
sim = sim(1+warmup:end);

%% check for missing values
% -999 is used to denote missing values in observed data, but this is later
% scaled by area. Therefore we check for all negative values, and ignore those.
idx = find(obs >= 0); 

%% Find the constant e
e = mean(obs(idx))/100;

%% Apply constant and transform flows
obs = log(obs+e);
sim = log(sim+e);

%% Calculate metric
top = sum((sim(idx) - obs(idx)).^2);
bot = sum((obs(idx) - mean(obs(idx))).^2);
val = 1 - (top/bot);
end

