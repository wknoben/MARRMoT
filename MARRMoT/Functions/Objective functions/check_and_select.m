function [sim, obs, idx] = check_and_select(sim, obs, idx)
% CHECK_AND_SELECT checks that sim and obs have the same number of
% elements, then filters them based on the values in idx AND on the steps 
% where obs >= 0.

% Copyright (C) 2021 Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
%
% Out:
% obs       - time series of observations       [mx1]
% sim       - time series of simulations        [mx1]
% idx       - vector of indices used to subset  [mx1]
%

%% 1. Check size of inputs:
%  make sure inputs are vertical and have the same size
obs = obs(:);
sim = sim(:);
if ~size(obs) == size(sim)
    error('Time series not of equal size.')
end

%% 2. Get indices where obs >= 0
% -999 is opten used to denote missing values in observed data. Therefore
% we check for all negative values, and ignore those.
idx_exists = find(obs >= 0);
 
%% 3. Update those if needed with user-input indices
% update default indices if idx is given
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
    end
end

%% 4. Filter to only selected indices
obs = obs(idx);
sim = sim(idx);