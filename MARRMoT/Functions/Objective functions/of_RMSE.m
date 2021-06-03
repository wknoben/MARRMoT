function [val,idx] = of_RMSE(obs,sim,idx)
% of_RMSE Calculates the Root Mean Squared Error of simulated streamflow.
% Ignores time steps with negative flow values.
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

%% Calculate metric
val = sqrt(mean((obs-sim).^2));
end

