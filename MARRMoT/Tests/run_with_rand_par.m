function output = run_with_rand_par(model)

% Load the data
load MARRMoT_example_data.mat

% Create a climatology data input structure. 
input_climatology.precip   = data_MARRMoT_examples.precipitation;                   % Daily data: P rate  [mm/d]
input_climatology.temp     = data_MARRMoT_examples.temperature;                     % Daily data: mean T  [degree C]
input_climatology.pet      = data_MARRMoT_examples.potential_evapotranspiration;    % Daily data: Ep rate [mm/d]
input_climatology.delta_t  = 1;                                                     % time step size of the inputs: 1 [d]

% Initiate model object and set random parameter set
m = feval(model);
m.theta = m.parRanges(:,1) + rand(m.numParams,1).*(m.parRanges(:,2)-m.parRanges(:,1));

% Create vector of initial store values
input_s0 = zeros(1, m.numStores);
 
% Create struct for solver options
input_solver_opts.resnorm_tolerance = 1e-10;                                       % Root-finding convergence tolerance
input_solver_opts.rerun_maxiter   = 6;                                           % Maximum number of re-runs

m.input_climate = input_climatology;
m.solver_opts   = input_solver_opts;
m.S0            = input_s0;

% Run the model
[output.ex,...
 output.in,...
 output.ss,...
 output.wb,...
 output.solver_info] = ...
            m.get_output();