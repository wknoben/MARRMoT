function [ odeFunc,strFlux ] = createOdeApprox_IE(Sold,delta_t,varargin)
%createOdeApprox_IE Returns function that approximates an ODE through
%Implicit Euler.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% In:
% Sold     - storage at time = t
% delta_t  - time step size
% varargin - number of anonymous functions that is equal to the number of
%            stores
%
% Out:
% odeFunc  - numerical approximation of ODE's by Explicit Euler. This
%            function needs to be generated every time step
% strFlux  - string that specifies which storages values the model function
%            should use to update flux values per time step. For Explicit
%            Euler this is 'S(1,2,...n)old'. This function needs to be
%            generated only once

% IE: 
% (y(t+1)-y(t)) / delta_t = dy/dt
%
% In case of hydrological model:
% (S(t+1)-S(t)) / delta_t = modelFunction(S(t+1))
% (S(t+1)-S(t)) / delta_t - modelFunction(S(t+1)) = 0

% Sections:
% 1. Create an ODE approximation
% 2. Create a function that specifies which storage values we need to
% update fluxes. Here that would be the storages at t-1, i.e. S_old

%% 1. ODE approximation
% Check the number of equations we're making
numEq = length(Sold);

% Create a string that specifies the number of stores
str_store       = sprintf('S%i,',1:numEq);
str_store(end)  = [];                                                       % remove last comma

% Create a string for each store
str_full = sprintf(...
    ['@(',str_store,') [',...                                               % anon function call
    repmat(...                                                              % repeat the numerical approximation for the number of stores
    ['(S%i - Sold(%i))/delta_t - varargin{%i}(',str_store,');'],1,numEq),...% insert the proper variables in sprintf
    ']'],... 
    repelem(1:numEq,3));                                                    % insert variables for store numbers

% Evaluate everything into a single output
odeFunc = eval(str_full);

% Example of 'str_full'
% @(S1,S2,S3,S4,S5) 
%       [(S1 - Sold(1))/delta_t - varargin{1}(S1,S2,S3,S4,S5);
%        (S2 - Sold(2))/delta_t - varargin{2}(S1,S2,S3,S4,S5);
%        (S3 - Sold(3))/delta_t - varargin{3}(S1,S2,S3,S4,S5);
%        (S4 - Sold(4))/delta_t - varargin{4}(S1,S2,S3,S4,S5);
%        (S5 - Sold(5))/delta_t - varargin{5}(S1,S2,S3,S4,S5);]

%% 2. Which store values to update fluxes?
% Each model function uses a variable called "tmp_sFlux" of size
% (1,number_of_stores) to update the fluxes for each time step. Here we
% define which store values ought to be used for this time stepping scheme.
% This requires that we define a string that is evaluated outside this
% function, but within the model function, because in the case of implicit
% time-stepping schemes the required value (S(t+1)) has not yet been found.
%
% E.g. for Explicit Euler, flux(t+1) should be calculated using store(t)
% E.g. for Implicit Euler, flux(t+1) should be calculated using store(t+1)

% We're using "numEq" stores, no need to calculate that again
% numEq = length(Sold);

% Check if this output is requested and carry on if so
if nargout == 2

    % Create input string
    str_in = sprintf('tmp_sNew(%i),',1:numEq);
    str_in(end) = [];

    % Create full string that contains a function that assigns the appropriate
    % variables to "tmp_sFlux".
    strFlux = sprintf(['tmp_sFlux = [',str_in,'];']);

end

% Example of 'strFlux'
% 'tmp_sFlux = [tmp_sNew(1),tmp_sNew(2),tmp_sNew(3)];'

end

