 function [ sNew, fval, stopflag, stopiter ] = ...
                                 rerunSolver( solverName,...
                                              solverOptions,...
                                              solveFun,...
                                              maxIter,...
                                              resnorm_tolerance,...
                                              initGuess,...
                                              oldVal,...
                                              varargin )
%rerunSolver Restarts a root-finding solver with different starting points.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Initialize iteration counter, sampling checker and find number of ODEs
iter      = 1;
resnorm   = resnorm_tolerance + 1;                                           % i.e. greater than the required accuracy
numStores = length(initGuess);    
stopflag  = 1;                                                               % normal function run

% Initialise vector of sNew and fval for each iteration, this way you can
% keep the best one, not the last one.
sNew_v    = zeros(numStores, maxIter);
fval_v    = inf(numStores,maxIter);
resnorm_v = inf(1, maxIter);

% Create the constant parts of the PROBLEM structure
problem.solver      = solverName;                                           % I.e. 'fsolve' or 'lsqnonlin'
problem.options     = solverOptions;                                        % option structure from 'optimoptions'
problem.objective   = solveFun;                                             % function to be solved
if size(varargin,2) > 0; problem.lb = cell2mat(varargin(1)); end            % only relevant for 'lsqnonlin' solver
if size(varargin,2) > 1; problem.ub = cell2mat(varargin(2)); end            % as above

% Start the re-sampling
% Re-sampling uses different starting points for the solver to see if
% solution accuracy improves. Starting points are alternated as follows:
% 1. location where the solver got stuck
% 2. storages at previous time step
% 3. minimum values
% 4. maximum values
% 5. randomized values close to solution of previous time steps

while resnorm > resnorm_tolerance

    % Select the starting points
    switch iter
        case 3                                                             % 1. Low values (store minima or zero)
            if isempty(problem.lb)
                problem.x0 = zeros(numStores,1);
            else
                problem.x0 = problem.lb;
            end
        case 4                                                             % 2. High values (store maxima or 2E4)
            if isempty(problem.ub)
                problem.x0 = 2*10^4.*ones(numStores,1);
            else
                problem.x0 = min(2*10^4.*ones(numStores,1), problem.ub);
            end
        case 1
            problem.x0 = initGuess(:);                                     % 3. Location where solver got stuck
        case 2
            problem.x0 = oldVal(:);                                        % 4. Stores at t-1
        otherwise
            problem.x0 = max(zeros(numStores,1),...
                             oldVal(:)+(rand(numStores,1)-0.5));           % 5. Randomized values close to starting location
    end
   
    % Re-run the solver
    [sNew_v(:,iter), fval_v(:,iter), stopflag] = feval(solverName, problem);
    
    resnorm_v(iter) = sum(fval_v(:,iter).^2);
    [resnorm,stopiter] = min(resnorm_v);
    fval = fval_v(:,stopiter);
    sNew = sNew_v(:,stopiter);
    
    % Increase the iteration counter
    iter = iter + 1;
    
    % Break out of the loop of iterations exceed the specified maximum
    if iter >= maxIter
        stopflag = 0;                                                          % function stopped due to iteration count
        break
    end
end


