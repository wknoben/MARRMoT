function [val,c,w] = of_inverse_KGE(obs,sim,varargin)
% of_inverse_KGE Calculates Kling-Gupta Efficiency of the inverse of
% simulated streamflow (Gupta et al, 2009), intended to capture low flow
% aspects better (Pushpalatha et al, 2012). Ignores time steps with -999 
% values.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% varargin  - optional weights of components    [3x1]
%
% Out:
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         [3x1]
% w         - weights    [wr,wa,wb]             [3x1]
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
% Decomposition of the mean squared error and NSE performance criteria: 
% Implications for improving hydrological modelling. Journal of Hydrology, 
% 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
%
% Pushpalatha, R., Perrin, C., Moine, N. Le, & Andréassian, V. (2012). A 
% review of efficiency criteria suitable for evaluating low-flow 
% simulations. Journal of Hydrology, 420–421, 171–182. 
% https://doi.org/10.1016/j.jhydrol.2011.11.055

%% check inputs and set defaults
if nargin < 2
    error('Not enugh input arguments')
elseif nargin > 4
    error('Too many inputs.')    
end

% defaults
w = [1,1,1];
warmup = 0; % time steps to ignore when calculating     
    
% update defaults weights if needed  
if nargin == 3 || nargin == 4
    if min(size(varargin{1})) == 1 && max(size(varargin{1})) == 3           % check weights variable for size
        w = varargin{1};                                                    % apply weights if size = correct
    else
        error('Weights should be a 3x1 or 1x3 vector.')                     % or throw error
    end
end   

% update default warmup period if needed
if nargin == 4
    if size(varargin{2}) == [1,1]
        warmup = varargin{2};
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

% check that inputs are column vectors ('corr()' breaks with rows)
% obs and sim should have the same orientation when we reach here
if size(sim,1) < size(sim,2)
    sim = sim';
    obs = obs';
end

%% Apply warmup period
obs = obs(1+warmup:end);
sim = sim(1+warmup:end);

%% check for missing values
% -999 is used to denote missing values in observed data, but this is later
% scaled by area. Therefore we check for all negative values, and ignore those.
idx = find(obs >= 0);                                                    

%% invert the time series and add a small constant to avoid issues with 0 flows
% Pushpalatha et al (2012) suggests to set e at 1/100th of the mean of the
% observed flow, which is what we'll follow here. The constant is added
% before transforming flows.

% Find the constant
e = mean(obs)/100;

% Apply the constant and transform flows
obs = 1./(obs+e);
sim = 1./(sim+e);

%% calculate components
c(1) = corr(obs(idx),sim(idx));                                             % r: linear correlation
c(2) = std(sim(idx))/std(obs(idx));                                         % alpha: ratio of standard deviations
c(3) = mean(sim(idx))/mean(obs(idx));                                       % beta: bias 

%% calculate value
val = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2);    % weighted KGE

end