function [val,c,w] = of_KGE(obs,sim,varargin)
% of_KGE Calculates Kling-Gupta Efficiency of simulated streamflow (Gupta
% et al, 2009). Ignores time steps with -999 values.
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

%% check inputs and set defaults
if nargin < 2
    error('Not enugh input arguments')
elseif nargin > 3
    error('Too many inputs.')    
    
elseif nargin == 2
    w = [1,1,1];                                                            % no weights specified, use defaults
    
elseif nargin == 3
    if size(varargin{1}) == [1,3] | size(varargin{1}) == [3,1]              % check weights variable for size
        w = varargin{1};                                                    % apply weights if size = correct
    else
        error('Weights should be a 3x1 or 1x3 vector.')                     % or throw error
    end
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
c(1) = corr(obs(idx),sim(idx));                                             % r: linear correlation
c(2) = std(sim(idx))/std(obs(idx));                                         % alpha: ratio of standard deviations
c(3) = mean(sim(idx))/mean(obs(idx));                                       % beta: bias 

%% calculate value
val = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2);    % weighted KGE

end