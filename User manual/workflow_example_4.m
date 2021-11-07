% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Contact:  l.trotter@unimelb.edu.au

% This example workflow contains an example application of calibration of a
% single models to a single catchment.
% It includes 7 steps:
%
% 1. Data preparation
% 2. Model choice and setup
% 3. Model solver settings and time-stepping scheme
% 4. Calibration settings
% 5. Model calibration
% 6. Evaluation of calibration results
% 7. Output vizualization

% NOTE: this example uses a custom function 'my_cmaes' to perform the
% optimisation, it is a wrapper around 'cmaes' to ensure inputs and outputs
% are consistent with other MATLAB's optimisation algorithms (e.g.
% 'fminsearch' or 'fminsearchbnd').
% While the wrapper is provided as part of MARRMoT, it requires the source 
% code to 'cmaes' to function, it is available at: 
% http://cma.gforge.inria.fr/cmaes.m
%
% The wrapper is necessary for the optimiser to function within the
% MARRMoT_model.calibrate method.
% Alternatively any model can be calibrated using any optimisation
% algorithm using the MARRMoT_model.calc_par_fitness method which returns
% the value of an objective function and can be used as input to an
% optimiser.

%% 1. Prepare data
% Load the data
load MARRMoT_example_data.mat

% Create a climatology data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_climatology.precip   = data_MARRMoT_examples.precipitation;                   % Daily data: P rate  [mm/d]
input_climatology.temp     = data_MARRMoT_examples.temperature;                     % Daily data: mean T  [degree C]
input_climatology.pet      = data_MARRMoT_examples.potential_evapotranspiration;    % Daily data: Ep rate [mm/d]
input_climatology.delta_t  = 1;                                                     % time step size of the inputs: 1 [d]

% Extract observed streamflow
Q_obs = data_MARRMoT_examples.streamflow;

%% 2. Define the model settings and create the model object

model     = 'm_37_hbv_15p_5s';                                             % Name of the model function (these can be found in Supporting Material 2)
m = feval(model);
parRanges = m.parRanges;                                                   % Parameter ranges
numParams = m.numParams;                                                   % Number of parameters
numStores = m.numStores;                                                   % Number of stores
input_s0  = zeros(numStores,1);                                            % Initial storages (see note in paragraph 5 on model warm-up)

%% 3. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance; users have reported differences in simulation accuracy (KGE scores) during calibration between Matlab and Octave for a given tolerance. In certain cases, Octave seems to require tigther tolerances to obtain the same KGE scores as Matlab does.
input_solver_opts.rerun_maxiter   = 6;                                           % Maximum number of re-runs

%% 4. Define calibration settings
% Settings for 'my_cmaes'
% the opts struct is made up of two fields, the names are hardcoded, so
% they cannot be changed:
%    .sigma0:     initial value of sigma
%    .cmaes_opts: struct of options for cmaes, see cmaes documentation
%                 or type cmaes to see list of options and default values

% starting sigma
optim_opts.insigma = .3*(parRanges(:,2) - parRanges(:,1));                 % starting sigma (this is default, could have left it blank)

% other options
optim_opts.LBounds  = parRanges(:,1);                                      % lower bounds of parameters
optim_opts.UBounds  = parRanges(:,2);                                      % upper bounds of parameters
optim_opts.PopSize  = 4 + floor(3*log(numParams));                         % population size (default)
optim_opts.TolX       = 1e-6 * min(optim_opts.insigma);                    % stopping criterion on changes to parameters 
optim_opts.TolFun     = 1e-4;                                              % stopping criterion on changes to fitness function
optim_opts.TolHistFun = 1e-5;                                              % stopping criterion on changes to fitness function
%optim_opts.MaxIter    = 5;                                                 % just do 5 iterations, to check if it works
optim_opts.SaveFilename      = 'wf_ex_4_cmaesvars.mat';                    % output file of cmaes variables
optim_opts.LogFilenamePrefix = 'wf_ex_4_';                                 % prefix for cmaes log-files

% Other useful options
    % change to true to run in parallel on a pool of CPUs (e.g. on a cluster)
optim_opts.EvalParallel = false;                                           
    % uncomment to restart r-times until the condition
% r = 2;
% optim_opts.Restarts  = r;
% optim_opts.RestartIf = 'fmin > -.8'; % OF is inverted, so this restarts
%                                        unless the OF (KGE here) > 0.8

% initial parameter set
par_ini = mean(parRanges,2);                                               % same as default value

% Choose the objective function
of_name      = 'of_KGE';                                                   % This function is provided as part of MARRMoT. See ./MARRMoT/Functions/Objective functions
weights      = [1,1,1];                                                    % Weights for the three KGE components

%
% Time periods for calibration.
% Indices of timestept to use for calibration, here we are using 1 year for
% warm-up, and 2 years for calibration, leaving the rest of the dataset for
% independent evaluation.
n = length(Q_obs);
warmup = 365;
cal_idx = (warmup+1):(warmup+365*2);
eval_idx = max(cal_idx):n;

%% 5. Calibrate the model
% MARRMoT model objects have a "calibrate" method that takes uses a chosen
% optimisation algorithm and objective function to optimise the parameter
% set. See MARRMoT_model class for details.

% first set up the model
m.input_climate = input_climatology;
%m.delta_t       = input_climatology.delta_t;         % unnecessary if input_climate already contains .delta_t
m.solver_opts   = input_solver_opts;
m.S0            = input_s0;

[par_opt,...                                                               % optimal parameter set
    of_cal,...                                                             % value of objective function at par_opt
    stopflag,...                                                           % flag indicating reason the algorithm stopped
    output] = ...                                                          % other info about parametrisation
              m.calibrate(...                                              % call the calibrate method of the model object
                          Q_obs,...                                        % observed streamflow
                          cal_idx,...                                      % timesteps to use for model calibration
                          'my_cmaes',...                                   % function to use for optimisation (must have same structure as fminsearch)
                          par_ini,...                                      % initial parameter estimates
                          optim_opts,...                                   % options to optim_fun
                          of_name,...                                      % name of objective function to use
                          1,...                                            % should the OF be inversed?
                          weights);                                        % additional arguments to of_name
                   
%% 6. Evaluate the calibrated parameters on unseen data
% Run the model with calibrated parameters, get only the streamflow
Q_sim = m.get_streamflow([],[],par_opt);  
             
% Compute evaluation performance
of_eval = feval(of_name,...                                                % Objective function name (here 'of_KGE')
                Q_obs,...                                                  % Observed flow during evaluation period
                Q_sim,...                                                  % Simulated flow during evaluation period, using calibrated parameters            
                eval_idx,...                                               % Indices of evaluation period
                weights);                                                  % KGE component weights
                
%% 7. Visualise the results
              
% Prepare a time vector
t = data_MARRMoT_examples.dates_as_datenum;

% Compare simulated and observed streamflow
figure('color','w'); 
    box on;
    hold all; 
    
    % Flows
    l(1) = plot(t,Q_obs,'k');
    l(2) = plot(t,Q_sim,'r');
    
    % Dividing line
    l(3) = plot([t(max(cal_idx)),t(max(cal_idx))],[0,170],'--b','linewidth',2);
    l(4) = plot([t(warmup),t(warmup)],[0,170],'--g','linewidth',2);
    
    % Legend & text       
    l = legend(l,'Q_{obs}','Q_{sim}','Cal // Eval', 'warmup // Cal','Location','northwest');
    title('Model calibration and evaluation results')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    
    txt_cal  = sprintf('Calibration period \nKGE = %.2f ',of_cal);
    txt_eval = sprintf('Evaluation period \nKGE = %.2f',of_eval);
    text(t(round(mean(cal_idx))),57,txt_cal,'fontsize',16,'HorizontalAlignment', 'center');
    text(t(round(mean(eval_idx))),57,txt_eval,'fontsize',16,'HorizontalAlignment', 'center');
    set(gca,'fontsize',16);
    
    % Other settings
    datetick;
    ylim([0,60])
    set(gca,'TickLength',[0.005,0.005])
