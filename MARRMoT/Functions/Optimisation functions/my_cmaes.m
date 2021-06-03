function [x,...                                                            % Solution, returned as a real vector or real array.
            fval,...                                                       % Objective function value at the solution, returned as a real number.
            exitflag,...                                                   % Reason cmaes stopped, returned as an integer.
            output] = ...                                                  % Information about the optimization process, returned as a structure
                      my_cmaes(fun,...                                     % Function to minimize, specified as a function handle or function name
                               x0,...                                      % Initial point, specified as a real vector or real array
                               opts)                                       % Optimization options, specified as a structure

% This is a wrapper around CMA-ES optimisation algorithm to make its inputs
% and outputs comparable with MATLAB's fminsearch

% This requires cmaes.m, which can be retrieved from:
% http://cma.gforge.inria.fr/

% Check whether cmaes is properly installed
if ~exist('cmaes','file')
    disp(['Optimizer ''cmaes'' not found. Visit the following ',...
          'link to download: ', 'http://cma.gforge.inria.fr/cmaes.m']);
    return
end

% Check that the opts struct contains the options needed for cmaes, set
% defaults otherwise
if ~isfield(opts, 'sigma0') || isempty(opts.sigma0)
    sigma0 = [];
    disp(['initial sigma not set: default value will be used',...
          'see cmaes.m for details']);
else
    sigma0 = opts.sigma0;
end

if ~isfield(opts, 'cmaes_opts') || isempty(opts.cmaes_opts)
    cmaes_opts = cmaes();
    disp(['CMA-ES options not set: default value will be used',...
          'see cmaes.m for details']);
else
    cmaes_opts = opts.cmaes_opts;
end

% if fun is a character, transform it into a handle
if ~ishandle(fun)
    fun = @(varargin) feval(fun, varargin{:});
end

% run CMA-ES
[x,...
    fval,...
    n_feval,...
    exitflag] = ...
                cmaes(fun,...
                      x0,...
                      sigma0,...
                      cmaes_opts);
output.funcCount = n_feval;
output.algorithm = 'Evolution Strategy with Covariance Matrix Adaptation (CMA-ES)';