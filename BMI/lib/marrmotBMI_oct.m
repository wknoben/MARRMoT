classdef marrmotBMI_oct < handle
    % MARRMoT m01 implemented as a BMI model
    
    properties
        % BMI required
        startTime           % start time
        endTime             % end time
        dt                  % time step size
        time_unit           % units of time
        time = 1;           % current time step
        lat                 % latitude
        lon                 % longitude
        
        % Indepdendent of model
        forcing = struct('precip',0,...             % time series of precipitation
                         'temp',0,...               % time series of temperature
                         'pet',0,...                % time series of potential evapotranspiration
                         'delta_t',0,...            % time step size in [days]
                         'time_unit','str');        % string with time units
        solver = struct('name','str',...            % string with solver function name
                        'resnorm_tolerance',0,...   % approximation accuracy
                        'resnorm_maxiter',0);       % max number of re-runs
        
            % Output storage for a single time-step
            output_externalFluxes = struct('Q',0,...    % Streamflow
                                           'Ea',0);     % Actual evap
            output_internalFluxes = struct('tmp',0);    % Internal fluxes, will be appended later
            output_modelStorages  = struct('S1',0);     % Model stores, can be appended later
            output_waterBalance   = 0;                  % Can be an optional argument to return a water balance check
                    
        % Dependent on the model
        model_name          % name of the model function file
        parameters          % parameter values from some source
        store_num           % number of stores
        path                % location of config file
        store_cur           % vector with current model states
    end
    
    
    %the following is the heart of the BMI class. Every BMI model needs to
    %implement these methods. For some methods a standard implementation is
    %provided (example: update_until). These can be overwritten by
    %implementing them in the model class file. See the AR1class.m for an
    %example. Methods that are not implemented will throw an error msg
    %declaring they are not implemented.
    %See the BMI reference for more details.
    methods
        
        % Unsure whether this is still needed - check
%         function obj=m01BMI()       
%         end
        
%% Model control functions

        % Initialize the model: this (1) loads data, (2) selects the model,
        % (3) sets all the settings needed to run the model; all from a
        % single _model-specific_ config file.
        function  initialize(obj,path)
            % Load config file
            load(path)                          % contains everything needed to run the model, plus BMI stuff

            % Assign values to object
            obj.model_name = model_name;      
            obj.startTime  = time_start;
            obj.endTime    = time_end;
            obj.dt         = forcing.delta_t_days;
            obj.time_unit  = forcing.time_unit;
            obj.forcing    = forcing;
            obj.solver     = solver;
            obj.parameters = parameters;
            obj.store_cur  = store_ini;         % this value gets updated on every update() call
            obj.store_num  = length(store_ini);
            obj.path       = path;
            obj.lat        = data_origin(1);
            obj.lon        = data_origin(2);
        end
        
        % Run model for 1 time step
        function update(obj)
           
           % Get the forcing for this time step
            input_forcing.precip       = obj.forcing.precip(obj.time);
            input_forcing.temp         = obj.forcing.temp(obj.time);
            input_forcing.pet          = obj.forcing.pet(obj.time);
            input_forcing.delta_t      = obj.forcing.delta_t_days;
                
           % Run model for 1 time step
            [output_ex,...                                                  % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
             output_in,...                                                  % Internal model fluxes
             output_ss,...                                                  % Internal storages
             ] = ...                                                        % Water balance check
                                feval(obj.model_name,...                    % Model function name
                                      input_forcing,...                     % Time series of climatic fluxes in simulation period
                                      obj.store_cur,...                     % Initial storages
                                      obj.parameters,...                    % Parameter values
                                      obj.solver);                          % Solver settings
            
            % Update variables
            obj.store_cur             = structfun(@mean,output_ss);         % ugly hack to get the current storages out of the structure, but works in octave
            obj.output_externalFluxes = output_ex;
            obj.output_internalFluxes = output_in;
            obj.output_modelStorages  = output_ss;
%            obj.output_waterBalance   = output_wb;
            obj.time                  = obj.time + obj.dt;
        end
        
        % Function to continously run the model
        function update_until(obj,until_time)
            while (obj.time <= until_time)
                obj.update;
            end
        end
        
        % Finalize (empty)
        function finalize()
        end
        
%% Variable getter and setter 

        % Get current values
        function output = get_value(obj,long_var_name)
            switch long_var_name
                case 'P'
                    output = obj.forcing.precip;
                case 'T'
                    output = obj.forcing.temp;
                case 'Ep'
                    output = obj.forcing.pet;
                case 'S(t)'
                    output = [obj.store_cur];
                case 'par'
                    output = [obj.parameters];
                case 'sol_resnorm_tolerance'
                    output = [obj.solver.resnorm_tolerance];
                case 'sol_resnorm_maxiter'
                    output = [obj.solver.resnorm_maxiter];
                case 'flux_out_Q'
                    output = [obj.output_externalFluxes.Q];
                case 'flux_out_Ea'
                    output = [obj.output_externalFluxes.Ea];
                % TODO: this model has no internal fluxes, to be implemented for other models
                % case 'flux_in_tmp'
                    % output = [obj.output_internalFluxes.tmp];
                case 'wb'
                    output = [obj.output_waterBalance];
                otherwise
                    error('unkown variable');
            end
        end
        
        function output = get_value_at_indices(obj,long_var_name, inds)
            if inds ~= [0,0]
                error(['Indices out of bounds.']);
            else
                output = get_value(obj,long_var_name);
            end
        end
        
        % Set inputs using handles
        function set_value(obj,long_var_name, src)
            switch long_var_name
                case 'P'
                    obj.forcing.precip = src;
                case 'T'
                    obj.forcing.temp = src;
                case 'Ep'
                    obj.forcing.pet = src;
                case 'S(t)'
                    obj.store_cur = src(1);
                case 'par'
                    obj.parameters = src(1);
                case 'sol_resnorm_tolerance'
                    obj.solver.resnorm_tolerance = src(1);
                case 'sol_resnorm_maxiter'
                    obj.solver.resnorm_maxiter = src(1);
                case 'flux_out_Q'
                    obj.output_externalFluxes.Q = src(1);
                case 'flux_out_Ea'
                    obj.output_externalFluxes.Ea = src(1);
                % TODO: this model has no internal fluxes, to be implemented for other models    
                % case 'tmp'
                %     obj.output_internalFluxes.tmp = src(1);
                case 'wb'
                    obj.output_waterBalance = src(1);
                otherwise
                    error('unkown variable');
            end
        end
        
        function set_value_at_indices(obj,long_var_name, inds, src)
            if inds ~= [0,0]
                error(['Indices out of bounds.']);
            else
                set_value(obj,long_var_name, src);
            end
        end
        
%% Model information functions        
        
        % Get input variable names
        function output = get_input_var_names(obj)
            % output={'Precipitation time series > P',... 
            %         'Temperature time series > T',...
            %         'Potential evapotranspiration time series > Ep',...
            %         'Storages (time = t) > S(t)',...
            %         'Model function name > mod',...
            %         'Parameters > par',...
            %         'Solver > sol'};
            output={
                'P', 'T', 'Ep', 'S(t)', 'par', 'sol_resnorm_tolerance',...
                'sol_resnorm_maxiter', 'flux_out_Q', 'flux_out_Ea',...
                'wb'
            }
        end
        
        % Get output variable names
        function output = get_output_var_names(obj)
            % output={'Fluxes leaving the model > flux_out',...
            %         'Internal fluxes > flux_in',...
            %         'Storages (time = t) > S(t)',...
            %         'Water balance check > wb'};
            output= {
                'P', 'T', 'Ep', 'S(t)', 'par', 'sol_resnorm_tolerance',...
                'sol_resnorm_maxiter', 'flux_out_Q', 'flux_out_Ea',...
                'wb'
            };
        end
        
        % Model components
        function output = get_component_name(obj)
            output=["MARRMoT ", obj.model_name, " ", obj.solver.name];
        end
        
%% Time functions        
        
        % Start time
        function output = get_start_time(obj)
            # datenum returns the date as a serial number since Jan 1 0000
            # output will return the date as a number of days since CUT - 1 Jan 1970 
            CUT_time = [1970 01 01 0 0 0];
            output = datenum(obj.startTime) - datenum(CUT_time);
            # output = datenum(obj.startTime);
        end
        
        % End time
        function output = get_end_time(obj)
            # output will return the date as a number of days since CUT - 1 Jan 1970 
            CUT_time = [1970 01 01 0 0 0];
            output = datenum(obj.endTime) - datenum(CUT_time);
        end
        
        % Time units
        function output = get_time_units(obj)
            timeformat_index = strncmp(obj.time_unit, 'days', 3);
            if timeformat_index ~= 0
              commandstring = "days since 1970-01-01 00:00:00.0 00:00";
              output = commandstring;
            else
              output = obj.time_unit;
            end
        end
        
        % Current time
        function output = get_current_time(obj)
            CUT_time = [1970 01 01 0 0 0];
            t_current = obj.time + datenum(obj.startTime)- datenum(CUT_time) - 1;
            output = t_current;
        end
        
        % Time step
        function output = get_time_step(obj)
            output=obj.dt;
        end

%% Model grid

        % Grid type
        function output = get_grid_type(~)
            output = 'uniform_rectilinear';         % Fixed for all MARRMoT models
        end

        % Grid rank
        function output = get_grid_rank(~)
            output = 2;                             % Fixed for all MARRMoT models
        end
        
        % Grid size
        function output = get_grid_size(~)
            output = 1;                             % Fixed for all MARRMoT models
        end
        
        % Grid shape
        function output = get_grid_shape(~)
            output = [1,1];                         % Fixed for all MARRMoT models
        end
        
        % Grid origin
        function output = get_grid_origin(obj)
            output = [obj.lat,obj.lon];             
        end
        
        % Grid spacing
        function output = get_grid_spacing(~)
            output = [0,0];                         % Fixed for all MARRMoT models
        end
        
        % Grid X
        function output = get_grid_x(~)
            output = 1;                              % Fixed for all MARRMoT models
        end
        
        % Grid Y
        function output = get_grid_y(~)
            output = 1;                              % Fixed for all MARRMoT models
        end
        
        % Grid Z
        function output = get_grid_z(~)
            output = 1;                              % Fixed for all MARRMoT models
        end

 
%% Variable information

        % Units
        function output = get_var_units(obj,long_var_name)
            switch long_var_name
                case 'P'
                    output = ['mm ',obj.time_unit];
                case 'T'
                    output = 'degree C';
                case 'Ep'
                    output = ['mm ',obj.time_unit];
                case 'S(t)'
                    output = 'mm';
                case 'par'
                    output = 'see model documentation';
                case 'sol_resnorm_tolerance'
                    output = '-';
                case 'sol_resnorm_maxiter'
                    output = '-';
                case 'flux_out_Q'
                    output = ['mm ',obj.time_unit];
                case 'flux_out_Ea'
                    output = ['mm ',obj.time_unit];
                case 'flux_in_tmp'
                    output = ['mm ',obj.time_unit];
                case 'wb'
                    output = 'mm';
                otherwise
                    error('unkown variable');
            end
        end
        
        % Type
        function output = get_var_type(obj,long_var_name)
            switch long_var_name
                case 'P'
                    tmp = obj.forcing.precip;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'T'
                    tmp = obj.forcing.temp;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'Ep'
                    tmp = obj.forcing.pet;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'S(t)'
                    tmp = obj.store_cur;
                    tmp = whos('tmp');
                    output = tmp.class;
                 case 'par'
                    tmp = obj.parameters;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'sol_resnorm_tolerance'
                    tmp = obj.solver.resnorm_tolerance;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'resnorm_maxiter'
                    tmp = obj.solver.resnorm_maxiter;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'flux_out_Q'
                    tmp = obj.output_externalFluxes.Q;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'flux_out_Ea'
                    tmp = obj.output_externalFluxes.Ea;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'flux_in_tmp'
                    tmp = obj.output_internalFluxes;
                    tmp = whos('tmp');
                    output = tmp.class;
                case 'wb'
                    tmp = obj.output_waterBalance;
                    tmp = whos('tmp');
                    output = tmp.class;
                otherwise
                    error('unkown variable');
            end
        end
        
        % Item size
        function output = get_var_itemsize(obj,long_var_name)
            switch long_var_name
                case 'P'
                    tmp = obj.forcing.precip;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'T'
                    tmp = obj.forcing.temp;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'Ep'
                    tmp = obj.forcing.pet;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'S(t)'
                    tmp = obj.store_cur;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'par'
                    tmp = obj.parameters;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'sol_resnorm_tolerance'
                    tmp = obj.solver.resnorm_tolerance;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'sol_resnorm_maxiter'
                    tmp = obj.solver.resnorm_maxiter;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'flux_out_Q'
                    tmp = obj.output_externalFluxes.Q;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'flux_out_Ea'
                    tmp = obj.output_externalFluxes.Ea;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'flux_in.tmp'
                    tmp = obj.output_internalFluxes.tmp;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                case 'wb'
                    tmp = obj.output_waterBalance;
                    tmp = whos('tmp');
                    output = tmp.bytes/(tmp.size(1)*tmp.size(2));
                otherwise
                    error('unkown variable');
            end
        end
        
        % Item nbytes
        function output = get_var_nbytes(obj,long_var_name)
            switch long_var_name
                case 'P'
                    tmp = obj.forcing.precip;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'T'
                    tmp = obj.forcing.temp;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'Ep'
                    tmp = obj.forcing.pet;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'S(t)'
                    tmp = obj.store_cur;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'par'
                    tmp = obj.parameters;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'sol_resnorm_tolerance'
                    tmp = obj.solver.resnorm_tolerance;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'sol_resnorm_maxiter'
                    tmp = obj.solver.resnorm_maxiter;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'flux_out_Q'
                    tmp = obj.output_externalFluxes.Q;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'flux_out_Ea'
                    tmp = obj.output_externalFluxes.Ea;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'flux_in_tmp'
                    tmp = obj.output_internalFluxes.tmp;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                case 'wb'
                    tmp = obj.output_waterBalance;
                    tmp = whos('tmp');
                    output = tmp.bytes;
                otherwise
                    error('unkown variable');
            end
        end

        % Grid
        function output = get_var_grid(obj,long_var_name)
            output = 1;                 % Grid size is 1 for all models, so arbitrary values
        end
        
        % Location
        function output = get_var_location(obj,long_var_name)
            output = 'face';            % Grid size is 1 for all models, so arbitrary values
        end
        
    end
end

