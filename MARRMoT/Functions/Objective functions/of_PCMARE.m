function [val,idx] = of_PCMARE(obs,sim,idx)
% of_MARE Calculates a version of the Mean Absolute Relative Error (MARE)
% of simulated streamflow as a percentage of the MARE of the mean 
% observed flow. Ignores time steps with negative flow values.  Adds a
% constant e of 1/100 of mean(obs) to avoid issues with zero flows 
% (Pushpalatha et al., 2012).

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
%
% Out:
% val       - objective function value          [1x1]
% idx       - vector of indices used for calculation

% Pushpalatha, R.; Perrin, C.; le Moine, N. and Andréassian V. (2012). "A
% review of efficiency criteria suitable for evaluating low-flow
% simulations". Journal of Hydrology. 420-421, 171-182. 
% doi:10.1016/j.jhydrol.2011.11.055

%% Check inputs and select timesteps
if nargin < 2
    error('Not enugh input arguments')    
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);                                             

%% Find the constant e
m = mean(obs);
e = m/100;

%% Apply constant and transform flows
obs = obs+e;
sim = sim+e;

%% Calculate metric
MARE = mean(abs((sim-obs)/obs));
MARE_mean = mean(abs((sim-m)/m));

val = MARE/MARE_mean;
end
