function [val, idx] = of_inverse_NSE(obs,sim,idx)
% of_inverse_NSE Calculates the Nash-Sutcliffe Efficiency (Nash & Sutcliffe, 1970) 
% of the log of simulated streamflow. Ignores time steps with negative flow
% values. Adds a constant e of 1/100 of mean(obs) to avoid issues with zero
% flows (Pushpalatha et al. 2012).

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
% idx       - indices used for the calculation
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

%% Check inputs and select timesteps
if nargin < 2
    error('Not enugh input arguments')    
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);

%% invert the time series and add a small constant to avoid issues with 0 flows
% Pushpalatha et al (2012) suggests to set e at 1/100th of the mean of the
% observed flow, which is what we'll follow here. The constant is added
% before transforming flows.

% Find the constant
e = mean(obs)/100;

% Apply the constant and transform flows
obs = 1./(obs+e);
sim = 1./(sim+e);

%% Calculate metric
top = sum((sim - obs).^2);
bot = sum((obs - mean(obs)).^2);
val = 1 - (top/bot);
end

