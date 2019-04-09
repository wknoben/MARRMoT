##% Script to run a MARRMoT model through the BMI interface
##
##%% 1. Create config file
##
##% Get some data
##load MARRMoT_example_data.mat
##
##% Catchment settings
##% -------------------------------------------------------------------------
##
##% 'data_origin' is a 2x1 array with [lat,lon]
##data_origin = [35.29,87.49];
##
##% 'forcing' is a structure with fields 'precip, 'temp', 'pet',
##% 'delta_t_days', 'time_unit'
##forcing.precip       = data_MARRMoT_examples.precipitation;
##forcing.temp         = data_MARRMoT_examples.temperature;
##forcing.pet          = data_MARRMoT_examples.potential_evapotranspiration;
##forcing.delta_t_days = 1;   % 1 [d]
##forcing.time_unit    = 'day';
##
##% 'time_start' & 'time_end' are a vector with starting/end date
##time_start = datevec(data_MARRMoT_examples.dates_as_datenum(1));
##time_end   = datevec(data_MARRMoT_examples.dates_as_datenum(end));
##
##% Model settings
##% -------------------------------------------------------------------------
##
##%'model_name' is a string with the name of the model function
##model_name = 'm_01_collie1_1p_1s';
##
##% 'parameters' is a vector of length n (number of parameters)
##parameters = 10; % Smax, [mm]
##
##% 'store_ini' is a vector of length m (number of stores)
##store_ini = 5; 
##
##% 'solver' is a structure with fields 'name', 'resnorm_tolerance',
##% 'resnorm_maxiter'
##solver.name = 'createOdeApprox_IE';
##solver.resnorm_tolerance = 0.1;
##solver.resnorm_maxiter = 6;
##
##% Save as a config file
##% -------------------------------------------------------------------------
##clear data_MARRMoT_examples
##save BMI_testcase_m01_BuffaloRiver_TN_USA
##
##    
##%% 2. Run the model through BMI
##clear all % start with a clean slate

addpath(genpath('../..'))

% Specify config path
filepath = '/Users/rwhut/Documents/eWaterCycle/repos/MARRMoT/BMI/Config/BMI_testcase_m01_BuffaloRiver_TN_USA';

% Create a model object
obj_marrmot_m01 = marrmotBMI_oct();

% Initialize the model
obj_marrmot_m01.initialize(filepath)

% Get the number of time steps
ts = datenum(obj_marrmot_m01.startTime);
te = datenum(obj_marrmot_m01.endTime);
tn = te-ts;

% Run a loop
qsim = NaN.*zeros(tn,1);

while obj_marrmot_m01.time <= tn
    
    % Display time step
    disp(obj_marrmot_m01.time)
    
    % Advance 1 time step
    obj_marrmot_m01.update();
    
    % Get simulated flow
    tmp = obj_marrmot_m01.get_value('flux_out'); % gives a structure
    qsim(obj_marrmot_m01.time-1) = tmp.Q;
    
end

%% 3. Do some visualisation

% get observations
load '/Users/rwhut/Documents/eWaterCycle/repos/MARRMoT/BMI/Config/MARRMoT_example_data.mat'

% Figure
figure('color','w');
hold on;
plot(data_MARRMoT_examples.streamflow);
plot(qsim)
title([obj_marrmot_m01.model_name, ' at ', num2str(obj_marrmot_m01.lat), ...
    ' lat, ',num2str(obj_marrmot_m01.lon), ' lon'],'interpreter','none')
xlabel('time [d]')
ylabel('discharge [mm/d]')



















