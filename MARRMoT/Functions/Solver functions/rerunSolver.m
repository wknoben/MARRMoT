 function [Snew, fval, stopflag, iter ] = ...
                                 rerunSolver( solverName,...
                                              solverOptions,...
                                              solveFun,...
                                              maxIter,...
                                              resnorm_tolerance,...
                                              initGuess,...
                                              oldVal,...
                                              varargin )
%rerunSolver Restarts a root-finding solver with different starting points.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Initialize iteration counter, sampling checker and find number of ODEs
iter      = 1;
resnorm   = resnorm_tolerance + 1;                                           % i.e. greater than the required accuracy
numStores = length(initGuess);    
stopflag  = 1;                                                               % normal function run

% % Initialise vector of sNew and fval for each iteration, this way you can
% % keep the best one, not the last one.
% Snew_v    = zeros(numStores, maxIter);
% fval_v    = inf(numStores,maxIter);
% resnorm_v = inf(1, maxIter);
% Snew = -1 * ones(numStores, 1);

% Create the constant parts of the PROBLEM structure
problem.solver      = solverName;                                           % I.e. 'fsolve' or 'lsqnonlin'
problem.options     = solverOptions;                                        % option structure from 'optimoptions'
problem.objective   = solveFun;                                             % function to be solved

if size(varargin,2) > 0; lb = cell2mat(varargin(1)); else; lb = zeros(numStores,1); end
if size(varargin,2) > 1; ub = cell2mat(varargin(2)); else; ub = inf(numStores,1); end
problem.lb = lb(:);
problem.ub = ub(:);

% Start the re-sampling
% Re-sampling uses different starting points for the solver to see if
% solution accuracy improves. Starting points are alternated as follows:
% 1. location where the solver got stuck
% 2. storages at previous time step
% 3. minimum values
% 4. maximum values
% 5. randomized values close to solution of previous time steps

while resnorm > resnorm_tolerance
    
    switch iter
        case 1
            problem.x0 = zeros(numStores,1);                                           % 1. Low values
        case 2
            problem.x0 = 2*10^4.*ones(numStores,1);                                    % 2. High values
        case 3
            problem.x0 = initGuess(:);                                                 % 3. Location where solver got stuck
        case 4
            problem.x0 = oldVal(:);                                                    % 4. Storages at t-1
        case 5
            if ~isinf(varargin{2})
                problem.x0 = ub(:);                                           % 5. Store maximum values (not always provided)
            else
                problem.x0 = max(zeros(numStores,1),...
                             oldVal(:)+(rand(numStores,1)-0.5));            % 6. Randomized values close to starting location
            end
        otherwise
            problem.x0 = max(zeros(numStores,1),...
                             oldVal(:)+(rand(numStores,1)-0.5));            % 7. Randomized values close to starting location
    end
   
    % Re-run the solver
    [Snew, resnorm, fval, stopflag] = feval(solverName, problem);
        
    % Increase the iteration counter
    iter = iter + 1;
    
    % Break out of the loop of iterations exceed the specified maximum
    if iter >= maxIter
        stopflag = 0;                                                          % function stopped due to iteration count
        break
    end
end


