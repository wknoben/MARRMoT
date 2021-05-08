function [val] = of_NRMSE(obs,sim,varargin)
% of_NRMSE Calculates the Normalized Root Mean Squared Error of simulated 
% streamflow, with normalization by the range of observations. Ignores time
% steps with negative flow values.
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


%% Calculate metric
norm = max(obs(idx))-min(obs(idx));
val = sqrt(mean((obs(idx)-sim(idx)).^2)) / norm;
end

