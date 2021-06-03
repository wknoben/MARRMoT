function [val,c,idx,w] = of_mean_hilo_root5_KGE(obs,sim,idx,w)
% of_mean_hilo_KGE Calculates the average Kling-Gupta Efficiency of the 
% normal and fifth root of simulated streamflow (Gupta et al, 2009).
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% In:
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
% w         - optional weights of components    {[3x1] [3x1]} {high low}
                                            %or [3x1] high = low
%
% Out:
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         {[3x1] [3x1]} {high low}
% idx       - optional vector of indices to use for calculation, can be
%               logical vector [nx1] or numeric vector [mx1], with m <= n
% w         - weights    [wr,wa,wb]             {[3x1] [3x1]} {high low}
%
% Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
% Decomposition of the mean squared error and NSE performance criteria: 
% Implications for improving hydrological modelling. Journal of Hydrology, 
% 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003


%% do some basic imput checking
if nargin < 3 || isempty(idx)
    idx = 1:numel(obs);
end

if nargin < 4 || isempty(w)
    w1 = [1 1 1];
    w2 = [1 1 1];
else
    if iscell(w) && size(w) == 2
        w1 = w{1}; w2 = w{2};
    elseif isarray(w)
        w1 = w; w2 = w;
    else
        error(['w should be either'...
               'a cell array of size 2, or'...
               'an array of size [1,3]']);
    end
end

%% call individual KGE functions (they have their own error checking)
[val1,c1,~  ,w1] = of_KGE(obs,sim,idx,w1);
[val2,c2,idx,w2] = of_inverse_KGE(obs,sim,idx,w2);

%% calculate value
val = 0.5*(val1+val2);      % weighted KGE
c   = {c1 c2};              % components from high and low KGE
w   = {w1 w2};              % weights from high and low KGE

end