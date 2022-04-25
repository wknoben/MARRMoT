function output = run_with_rand_par(model)

% Load the data
load MARRMoT_example_data.mat

% Create a climatology data input structure. 
input_climatology.precip   = data_MARRMoT_examples.precipitation;                   % Daily data: P rate  [mm/d]
input_climatology.temp     = data_MARRMoT_examples.temperature;                     % Daily data: mean T  [degree C]
input_climatology.pet      = data_MARRMoT_examples.potential_evapotranspiration;    % Daily data: Ep rate [mm/d]

delta_t  = 1;                                                              % time step size of the inputs: 1 [d]

% Initiate model object and set random parameter set
m = feval(model, delta_t);
m.theta = m.parRanges(:,1) + rand(m.numParams,1).*(m.parRanges(:,2)-m.parRanges(:,1));

% Create vector of initial store values
input_s0 = zeros(1, m.numStores);
 
% Create struct for solver options
solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
solver_opts.resnorm_maxiter   = 6;                                         % Maximum number of re-runs
solver_opts.NewtonRaphson = optimset('MaxIter', m.numStores * 10);
solver_opts.fsolve = optimoptions('fsolve',...
                                  'Display','none',...                     % Disable display settings
                                  'JacobPattern', m.JacobPattern);
solver_opts.lsqnonlin = optimoptions('lsqnonlin',...                       % lsqnonlin settings for cases where fsolve fails
                                     'Display','none',...
                                     'JacobPattern',m.JacobPattern,...
                                      'MaxFunEvals',1000);
                                  
% Run the model
[output.ex,...
 output.in,...
 output.ss] = ...
            m.get_output(input_climatology,...
                         input_s0,...
                         solver_opts);