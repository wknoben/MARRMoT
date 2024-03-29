% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Contact:  l.trotter@unimelb.edu.au

% This example workflow contains an example application of a single model 
% to a single catchment.
% It includes 5 steps:
%
% 1. Data preparation
% 2. Model choice and setup
% 3. Model solver settings
% 4. Model generation and set-up
% 5. Model runs
% 6. Output vizualization

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

%% 2. Define the model settings
% NOTE: this example assumes that parameter values for this combination of
% model and catchment are known. 

% Model name 
% NOTE: these can be found in the Model Descriptions.
model       = 'm_29_hymod_5p_5s';                     

% Initial storage values
% NOTE: see the model function for the order in which stores are given. For
% HyMOD, this is on lines 86-91.

input_s0       = [ 1;                                                       % Initial soil moisture storage [mm]
                   7;                                                       % Initial fast flow 1 storage [mm]
                   3;                                                       % Initial fast flow 2 storage [mm] 
                   8;                                                       % Initial fast flow 3 storage [mm] 
                  22];                                                      % Initial slow flow storage [mm]

%% 3. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
input_solver_opts.resnorm_maxiter   = 6;                                           % Maximum number of re-runs
% these are the same settings that run by default if no settings are given

%% 4. Create a model object
% Create a model object
m = feval(model);

% Set up the model
m.input_climate = input_climatology;
%m.delta_t       = input_climatology.delta_t;         % unnecessary if input_climate already contains .delta_t
m.solver_opts   = input_solver_opts;
m.S0            = input_s0;


% Extract parameter ranges
model_range = m.parRanges;

%% 5. Run the model and extract all outputs
% Define the requested number of samples
numSample = 10;
numPar    = m.numParams;

% Prepare a storage array
results_mc_sampling = cell(numSample+1,5);

% Add a header
results_mc_sampling{1,1} = 'parameter_values';
results_mc_sampling{1,2} = 'output_ex';
results_mc_sampling{1,3} = 'output_in';
results_mc_sampling{1,4} = 'output_ss';
results_mc_sampling{1,5} = 'output_wb';

% Start the sampling
for i = 1:numSample

    % Sample a parameter set from the range
    input_theta = model_range(:,1)+rand(numPar,1).*(model_range(:,2)-model_range(:,1));
    
    % Run the model with the given parameter set
    [output_ex,...                                                             % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
     output_in,...                                                             % Internal model fluxes
     output_ss ,...                                                            % Internal storages
     output_waterbalance] = ...                                                % Water balance check              
                   m.get_output([],[],...                                      % Model method to run and return all outputs
                                input_theta);
    % Save the results
    results_mc_sampling{1+i,1} = input_theta;
    results_mc_sampling{1+i,2} = output_ex;
    results_mc_sampling{1+i,3} = output_in;
    results_mc_sampling{1+i,4} = output_ss;
    results_mc_sampling{1+i,5} = output_waterbalance;
    
    % Display a separation line
    disp(' ')
    disp('-------------------------------------')
    disp(' ')
    
end    
    
%% 5. Analyze the outputs                   
% Prepare a time vector
t = data_MARRMoT_examples.dates_as_datenum;

% Compare simulated and observed streamflow
figure('color','w'); 
    box on;
    hold on; 
    
    for i = 1:numSample
        h(1) = plot(t,results_mc_sampling{1+i,2}.Q,'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,data_MARRMoT_examples.streamflow,'r:');
    
    legend(h,'Simulated','Observed')
    title('Monte Carlo sampling results')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);

clear h i t 
