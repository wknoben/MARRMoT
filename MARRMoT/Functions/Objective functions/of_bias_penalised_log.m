function [val, idx] = of_bias_penalised_log(obs,sim, idx, of_name,varargin)
% of_bias_penalised_log applies a bias penalisation to an objective
% function as suggested by Viney et al. (2009). Ignore timesteps with
% negative values.
%
% Copyright (C) 2021 L. Trotter
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Given the efficiency of te OF specified (E_of):
%       E = E_of - 5 * abs(ln(1 + B)) ^ 2.5
% where B = mean(sim - obs) / mean(obs) is the bias between obs and sim
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
% of_name   - string name of the objective function
% varargin  - additional arguments to of_name
%
% Out:
% val       - objective function value          [1x1]
% idx       - indices used for the calculation
%
% Viney, N.R., Perraud, J., Vaze, J., Chiew, F.H.S., Post, D.A., Yang, A., 
% 2009. The usefulness of bias constraints in model calibration for 
% regionalisation to ungauged catchments, in: Proceedings of the 18 Th
% World IMACS / MODSIM Congress. Cairns, Australia.

%% check inputs
if nargin < 2 || isempty(of_name)
    error('Not enugh input arguments')
end

if nargin < 3; idx = []; end
[sim, obs, idx] = check_and_select(sim, obs, idx);

%% calculate components
E_of = feval(of_name, obs, sim, [], varargin{:});
B = mean(sim - obs)/mean(obs);

%% calculate value

val = E_of - 5 * abs(log(1 + B)) ^ 2.5;