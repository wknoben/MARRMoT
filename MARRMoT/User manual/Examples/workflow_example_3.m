% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Contact:  l.trotter@unimelb.edu.au

% This example workflow  contains an example application of multiple models
% to a single catchment.
% It includes 6 steps:
%
% 1. Data preparation
% 2. Model choice and setup
% 3. Model solver settings
% For each model in list
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
% NOTE: this example assumes that the model parameters for each model will 
% be sampled as part of the investigation. See lines 69-77 in this script.

% Model name 
% NOTE: these can be found in the Model Descriptions.
model_list  = {'m_29_hymod_5p_5s',...
               'm_01_collie1_1p_1s',...
               'm_27_tank_12p_4s'};                     

%% 3. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
input_solver_opts.rerun_maxiter   = 6;                                           % Maximum number of re-runs
% these are the same settings that run by default if no settings are given


%% Create and run all model objects
% Prepare a storage array
results_sampling = cell(length(model_list)+1,6);

% Add a header
results_sampling{1,1} = 'model name';
results_sampling{1,2} = 'parameter_values';
results_sampling{1,3} = 'output_ex';
results_sampling{1,4} = 'output_in';
results_sampling{1,5} = 'output_ss';
results_sampling{1,6} = 'output_wb';

% Start the sampling
for i = 1:length(model_list)
    
    % Define the current model
    model = model_list{i};
    
    % Display progress
    % ---------------------------------------------------------------------
    disp(['Now starting model ',model,'.'])
    disp(' ')
   
    %% 4. Create each model object    
    % Create a model object
    m = feval(model);

    % Extract parameter ranges
    model_range = m.parRanges;

    %% 5. Run the model and extract all outputs
    
    % Find number of stores and parameters
    numPar      = m.numParams;
    numStore    = m.numStores;
    
    % Sample a parameter set from the range
    input_theta = model_range(:,1)+rand(numPar,1).*(model_range(:,2)-model_range(:,1));
    
    % Set the inital storages
    input_s0    = zeros(numStore,1);
    
    
    % Run the model
    [output_ex,...                                                             % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
     output_in,...                                                             % Internal model fluxes
     output_ss ,...                                                            % Internal storages
     output_waterbalance] = ...                                                % Water balance check              
                   m.get_output(...                                            % Model method to run and return all outputs
                                input_climatology,...                          % Time series of climatic fluxes in simulation period
                                input_s0,...                                   % Initial storages
                                input_theta,...                                % Parameter set.
                                input_solver_opts);                            % Options for numerical solving of ODEs
    
    % Save the results
    % ---------------------------------------------------------------------
    results_sampling{1+i,1} = model;
    results_sampling{1+i,2} = input_theta;
    results_sampling{1+i,3} = output_ex;
    results_sampling{1+i,4} = output_in;
    results_sampling{1+i,5} = output_ss;
    results_sampling{1+i,6} = output_waterbalance;
    
    % Display a separation line
    disp(' ')
    disp('-------------------------------------')
    disp(' ')
    
end    
    
%% 6. Analyze the outputs                   
% Prepare a time vector
t = data_MARRMoT_examples.dates_as_datenum;

% Compare simulated and observed streamflow
figure('color','w'); 
    box on;
    hold all; 
    
    h(1) = plot(t,data_MARRMoT_examples.streamflow,'k','linewidth',2);
    for i = 1:length(model_list)
        h(1+i) = plot(t,results_sampling{1+i,3}.Q);
    end
        
    legend(h,['Observed',model_list], 'interpreter','none');
    title('Model sampling results')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    ylim([0,70])
    set(gca,'fontsize',16);

clear h i lh t 
    











