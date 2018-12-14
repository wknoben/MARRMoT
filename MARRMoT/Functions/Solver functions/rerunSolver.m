 function [ tmp_sNew, resnorm, flag ] = rerunSolver( solverName, solverOptions, solveFun, maxIter, minX, initGuess, oldVal, varargin )
%rerunSolver Restarts a root-finding solver with different starting points.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Initialize iteration counter, sampling checker and find number of ODEs
iter    = 1;
resnorm = minX + 1;                                                         % i.e. greater than the required accuracy
numODE  = length(initGuess);    
flag    = 0;                                                                % normal function run

% Create the constant parts of the PROBLEM structure
problem.solver      = solverName;                                           % I.e. 'fsolve' or 'lsqnonlin'
problem.options     = solverOptions;                                        % option structure from 'optimoptions'
problem.objective   = solveFun;                                             % function to be solved
if size(varargin,2) > 0; problem.lb = cell2mat(varargin(1)); end            % only relevant for 'lsqnonlin' solver
if size(varargin,2) > 1; problem.ub = cell2mat(varargin(2)); end            % as above


% Start the re-sampling
% Re-sampling uses different starting points for the solver to see if
% solution accuracy improves. Starting points are alternated as follows:
% 1. minimum values
% 2. maximum values
% 3. location where the solver got stuck
% 4. randomized values close to solution of previous time steps

while resnorm > minX

    % Select the starting points
    switch iter
        case 1
            problem.x0 = zeros(numODE,1);                                           % 1. Low values
        case 2
            problem.x0 = 2*10^4.*ones(numODE,1);                                    % 2. High values
        case 3
            problem.x0 = initGuess;                                                 % 3. Location where solver got stuck
        case 4
            problem.x0 = oldVal;                                                    % 4. Storages at t-1
        case 5
            if ~isempty(varargin{2})
                problem.x0 = varargin{2};                                           % 5. Store maximum values (not always provided)
            else
                problem.x0 = max(zeros(1,numODE),oldVal+(rand(1,numODE)-0.5));      % 6. Randomized values close to starting location
            end
        otherwise
            problem.x0 = max(zeros(1,numODE),oldVal+(rand(1,numODE)-0.5));          % 6. Randomized values close to starting location
    end
   
    % Re-run the solver
    [tmp_sNew,resnorm] = feval(solverName, problem);
    
    % Increase the iteration counter
    iter = iter + 1;
    
    % Break out of the loop of iterations exceed the specified maximum
    if iter >= maxIter
        flag = -1;                                                          % function stopped due to iteration count
        break
    end
end


