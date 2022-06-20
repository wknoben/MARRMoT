function [x, F, exitflag] = NewtonRaphson_octave(fun, x0, options)
% NEWTONRAPHSON Solve set o non-linear equations using Newton-Raphson method.
% this version works in Octave

% Copyright (C) 2021 Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This code is modified from the original found at:
% https://github.com/mikofski/NewtonRaphson

% Copyright (c) 2021, Mark Mikofski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% [X, RESNORM, F, EXITFLAG] = NEWTONRAPHSON(FUN, X0, OPTIONS)
% FUN is a function handle that returns a vector of residuals equations, F,
% and takes a vector, x, as its only argument. When the equations are
% solved by x, then F(x) == zeros(size(F(:), 1)).
%
% Optionally FUN may return the Jacobian, Jij = dFi/dxj, as an additional
% output. The Jacobian must have the same number of rows as F and the same
% number of columns as x. The columns of the Jacobians correspond to d/dxj and
% the rows correspond to dFi/d.
%
%   EG:  J23 = dF2/dx3 is the 2nd row ad 3rd column.
%
% If FUN only returns one output, then J is estimated using a center
% difference approximation,
%
%   Jij = dFi/dxj = (Fi(xj + dx) - Fi(xj - dx))/2/dx.
%
% NOTE: If the Jacobian is not square the system is either over or under
% constrained.
%
% X0 is a vector of initial guesses.
%
% OPTIONS is a structure of solver options created using OPTIMSET.
% EG: options = optimset('TolX', 0.001).
%
% The following options can be set:
% * OPTIONS.TOLFUN is the maximum tolerance of the norm of the residuals.
%   [1e-6]
% * OPTIONS.TOLX is the minimum tolerance of the relative maximum stepsize.
%   [1e-6]
% * OPTIONS.MAXITER is the maximum number of iterations before giving up.
%   [100]
%
% X is the solution that solves the set of equations within the given tolerance.
% RESNORM is norm(F) and F is F(X). EXITFLAG is an integer that corresponds to
% the output conditions, OUTPUT is a structure containing the number of
% iterations, the final stepsize and exitflag message and JACOB is the J(X).
%
% See also OPTIMSET, OPTIMGET, FMINSEARCH, FZERO, FMINBND, FSOLVE, LSQNONLIN
%
% References:
% * http://en.wikipedia.org/wiki/Newton's_method
% * http://en.wikipedia.org/wiki/Newton's_method_in_optimization
% * 9.7 Globally Convergent Methods for Nonlinear Systems of Equations 383,
%   Numerical Recipes in C, Second Edition (1992),
%   http://www.nrbook.com/a/bookcpdf.php
% Version 0.5
% * allow sparse matrices, replace cond() with condest()
% * check if Jstar has NaN or Inf, return NaN or Inf for cond() and return
%   exitflag: -1, matrix is singular.
% * fix bug: max iteration detection and exitflag reporting typos
% Version 0.4
% * allow lsq curve fitting type problems, IE non-square matrices
% * exit if J is singular or if dx is NaN or Inf
% Version 0.3
% * Display RCOND each step.
% * Replace nargout checking in funwrapper with ducktypin.
% * Remove Ftyp and F scaling b/c F(typx)->0 & F/Ftyp->Inf!
% * User Numerical Recipies minimum Newton step, backtracking line search
%   with alpha = 1e-4, min_lambda = 0.1 and max_lambda = 0.5.
% * Output messages, exitflag and min relative step.
% Version 0.2
% * Remove `options.FinDiffRelStep` and `options.TypicalX` since not in MATLAB.
% * Set `dx = eps^(1/3)` in `jacobian` function.
% * Remove `options` argument from `funwrapper` & `jacobian` functions
%   since no longer needed.
% * Set typx = x0; typx(x0==0) = 1; % use initial guess as typx, if 0 use 1.
% * Replace `feval` with `evalf` since `feval` is builtin.
%% initialize
% There are no argument checks!
x0 = x0(:); % needs to be a column vector
% set default options

defaultopt = struct('TolX', 1e-12, 'TolFun', 1e-6, 'MaxIter', 1000);%, 'Display', 'off');

%% get options
TOLX = optimget(options, 'TolX', defaultopt.TolX);
TOLFUN = optimget(options, 'TolFun', defaultopt.TolFun);
MAXITER = optimget(options, 'MaxIter', defaultopt.MaxIter);

%TYPX = max(abs(x0), 1); % x scaling value, remove zeros
ALPHA = 1e-4; % criteria for decrease
MIN_LAMBDA = 0.1; % min lambda
MAX_LAMBDA = 0.5; % max lambda
%% set scaling values
% TODO: let user set weights
%weight = ones(numel(fun(x0)),1);
%J0 = weight*(1./TYPX'); % Jacobian scaling matrix
%% check initial guess
x = x0; % initial guess
F = fun(x); % evaluate initial guess
nf = length(F);
J = jacobian(fun, x, nf, F);
%Jstar = J./J0; % scale Jacobian
%if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
if any(isnan(J(:))) || any(isinf(J(:)))
    exitflag = -1; % matrix may be singular
else
    exitflag = 1; % normal exit
end
resnorm = norm(F, Inf); % calculate norm of the residuals
resnorm0 = 100*resnorm;dx = zeros(size(x0)); % dummy values
%% solver
Niter = 0; % start counter
lambda = 1; % backtracking
while (resnorm>TOLFUN || lambda<1) && exitflag>=0 && Niter<=MAXITER
    if lambda==1
        %% Newton-Raphson solver
        Niter = Niter+1; % increment counter
        % update Jacobian, only if necessary
        if resnorm/resnorm0 > .2
            J = jacobian(fun, x, nf, F);
%            Jstar = J./J0; % scale Jacobian
%            if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
            if any(isnan(J(:))) || any(isinf(J(:)))
                exitflag = -1; % matrix may be singular
                break
            end
        end
%        dx_star = -Jstar\F; % calculate Newton step
%        % NOTE: use isnan(f) || isinf(f) instead of STPMAX
%        dx = dx_star.*TYPX; % rescale x
        if rcond(J) <= eps; dx = pinv(J) * -F; else dx = -J\F; end
        g = F'*J;%star; % gradient of resnorm
        slope = g*dx;%_star; % slope of gradient
        fold = F'*F; % objective function
        xold = x; % initial value
        lambda_min = TOLX/max(abs(dx)./max(abs(xold), 1));
    end
    if lambda<lambda_min
        exitflag = 2; % x is too close to XOLD
        break
    elseif any(isnan(dx)) || any(isinf(dx))
        exitflag = -1; % matrix may be singular
        break
    end
    x = xold+dx*lambda; % next guess
    F = fun(x); % evaluate this guess
    f = F'*F; % new objective function
    %% check for convergence
    lambda1 = lambda; % save previous lambda
    if f>fold+ALPHA*lambda*slope
        if lambda==1
            lambda = -slope/2/(f-fold-slope); % calculate lambda
        else
            A = 1/(lambda1 - lambda2);
            B = [1/lambda1^2,-1/lambda2^2;-lambda2/lambda1^2,lambda1/lambda2^2];
            C = [f-fold-lambda1*slope;f2-fold-lambda2*slope];
            coeff = num2cell(A*B*C);
            [a,b] = coeff{:};
            if a==0
                lambda = -slope/2/b;
            else
                discriminant = b^2 - 3*a*slope;
                if discriminant<0
                    lambda = MAX_LAMBDA*lambda1;
                elseif b<=0
                    lambda = (-b+sqrt(discriminant))/3/a;
                else
                    lambda = -slope/(b+sqrt(discriminant));
                end
            end
            lambda = min(lambda,MAX_LAMBDA*lambda1); % minimum step length
        end
    elseif isnan(f) || isinf(f)
        % limit undefined evaluation or overflow
        lambda = MAX_LAMBDA*lambda1;
    else
        lambda = 1; % fraction of Newton step
    end
    if lambda<1
        lambda2 = lambda1;f2 = f; % save 2nd most previous value
        lambda = max(lambda,MIN_LAMBDA*lambda1); % minimum step length
        continue
    end
    resnorm0 = resnorm; % old resnorm
    resnorm = norm(F, Inf); % calculate new resnorm
end
%% output
% output.iterations = Niter; % final number of iterations
% output.stepsize = dx; % final stepsize
% output.lambda = lambda; % final lambda
if Niter>=MAXITER
    exitflag = 0;
%     output.message = 'Number of iterations exceeded OPTIONS.MAXITER.';
% elseif exitflag==2
%     output.message = 'May have converged, but X is too close to XOLD.';
% elseif exitflag==-1
%     output.message = 'Matrix may be singular. Step was NaN or Inf.';
% else
%     output.message = 'Normal exit.';
end
% jacob = J;
end
function J = jacobian(fun, x, nf, funx)
% estimate J
dx = eps^(1/3); % finite difference delta
nx = numel(x); % degrees of freedom
if nargin <4
    funx = fun(x); 
    if nargin <3
        nf = numel(funx);
    end % number of functions
end
J = zeros(nf,nx); % matrix of zeros
for n = 1:nx
    % create a vector of deltas, change delta_n by dx
    delta = J(:, n); delta(n) = dx;
    dF = fun(x+delta)-funx;%-delta); % delta F
    J(:, n) = dF(:)/dx;%/2; % derivatives dF/d_n
end
end