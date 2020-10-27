function [val] = of_bias_penalised_log(obs,sim,of_name,varargin)

% of_bias_penalised_log applies a bias penalisation to an objective
% function as suggested by Viney et al. (2009). Ignore timesteps with
% negative values.
%
% Given the efficiency of te OF specified (E_of):
%       E = E_of - 5 * abs(ln(1 + B)) ^ 2.5
% where B = mean(sim - obs) / mean(obs) is the bias between obs and sim
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% of_name   - string name of the objective function
% varargin  - additional arguments to of_name
%
% Out:
% val       - objective function value          [1x1]

%% check inputs
if nargin < 3
    error('Not enugh input arguments')
end

% check time series size and rotate one if needed
if checkTimeseriesSize(obs,sim) == 0
    error('Time series not of equal size.')
    
elseif checkTimeseriesSize(obs,sim) == 2
    sim = sim';                                                             % 2 indicates that obs and sim are the same size but have different orientations
end

%% check for missing values
% -999 is used to denote missing values in observed data, but this is later
% scaled by area. Therefore we check for all negative values, and ignore those.
idx = find(obs >= 0); 

%% calculate components
E_of = feval(of_name, obs, sim, varargin{:});
B = mean(sim(idx) - obs(idx))/mean(obs(idx));                                       % beta: bias 

%% calculate value

val = E_of - 5 * abs(log(1 + B)) ^ 2.5;