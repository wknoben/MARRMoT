classdef MARRMoT_model < handle
    % Class for MARRMoT models
    properties
        % static attributes, set for each models in the model definition
        numStores         % number of model stores
        numFluxes         % number of model fluxes
        numParams         % number of model parameters
        parRanges         % default parameter ranges
        JacobPattern      % pattern of the Jacobian matrix of model store ODEs
        StoreNames        % Names for the stores
        FluxNames         % Names for the fluxes
        FluxGroups        % Grouping of fluxes (useful for water balance and output)
        StoreSigns        % Signs to give to stores (-1 is a deficit store), assumes all 1 if not given
        % attributes set at the beginning of the simulation
            % directly by the user
        theta             % Set of parameters
        delta_t           % time step
        S0                % initial store values
        input_climate     % vector of input climate
        solver_opts       % options for numerical solving of ODEs
            % automatically, based on parameter set
        store_min         % store minimum values
        store_max         % store maximum values
        % attributes created and updated automatically throughout a
        % simulation
        t                 % current timestep
        fluxes            % vector of all fluxes
        stores            % vector of all stores
        uhs               % unit hydrographs
        fluxes_stf        % still-to-flow fluxes
        solver_data       % step-by-step info of solver used and residuals
        status            % 0 = model created, 1 = simulation ended
    end
    methods
        
        % Set methods with checks on inputs for attributes set by the user:
        function [] = set.delta_t(obj, value)
            if numel(value) == 1
                obj.delta_t = value;
                obj.reset();
            else
                error('delta_t must be a scalar')
            end
        end
        function [] = set.theta(obj, value)
            if numel(value) == obj.numParams
                obj.theta = value(:);
                obj.reset();
            else
                error(['theta must have ' int2str(obj.numParams) ' elements'])
            end
        end
        function [] = set.input_climate(obj, value)
            if isstruct(value)
                if isfield(value, 'delta_t')
                    obj.delta_t = value.delta_t;
                elseif isempty(obj.delta_t)
                    error(['delta_t is not in input climate struct: '...
                           'add it with obj.delta_t = delta_t'])
                end
                if isfield(value, {'precip' 'pet' 'temp'})
                    P = value.precip/obj.delta_t;
                    Ea = value.pet/obj.delta_t;
                    T = value.temp/obj.delta_t;
                    obj.input_climate = [P(:) Ea(:) T(:)];
                    obj.reset();
                else
                    error(['Input climate struct must contain fields: '...
                           'precip, pet, temp']);
                end
            elseif isnumeric(value)
                if size(value,2)
                    obj.input_climate = value;
                    obj.reset();
                else
                    error(['Input climate must have 3 columns: '...
                           'precip, pet, temp']);
                end
            else
                error(['Input climate must either be a struct or '...
                       'a numeric array of 3 columns']);
            end
        end
        function [] = set.S0(obj, value)
            if numel(value) == obj.numStores
                obj.S0 = value(:);
                obj.reset();
            else
                error(['S0 must have ' int2str(obj.numStores) ' elements'])
            end
        end
        function [] = set.solver_opts(obj, value)
            % add options to default ones
            obj.solver_opts = obj.add_to_def_opts(value);
            obj.reset();
        end
        
        % INIT_ runs before each model run to initialise store limits,
        % auxiliary parameters etc. it calls INIT which is model specific
        function obj = init_(obj)
            % min and max of stores
            obj.store_min = zeros(obj.numStores,1);
            obj.store_max = inf(obj.numStores,1);
            
            % empty vectors of fluxes and stores
            t_end = size(obj.input_climate, 1);
            obj.stores = zeros(t_end, obj.numStores);
            obj.fluxes = zeros(t_end, obj.numFluxes);
            
            % empty struct with the solver data
            obj.solver_data.resnorm   = zeros(t_end,1);
            obj.solver_data.solver = strings(t_end,1);
            obj.solver_data.iter   = zeros(t_end,1);
            
            % model specific initialisation
            obj.init();
        end
        
        % RESET is called any time that a user-specified input is changed
        % (t, delta_t, input_climate, S0, solver_options) and resets any
        % previous simulation ran on the object.
        % This is to prevent human error in analysing results.
        function obj = reset(obj)
        	obj.t = [];                 % current timestep
            obj.fluxes = [];            % vector of all fluxes
            obj.stores = [];            % vector of all stores
            obj.uhs = [];               % unit hydrographs
            obj.fluxes_stf = [];        % still-to-flow fluxes
            obj.solver_data = [];       % step-by-step info of solver used and residuals
            obj.status = 0;             % 0 = model created, 1 = simulation ended
        end
        
        % ODE approximation with Implicit Euler time-stepping scheme
        function err = solve_fun_IE(obj, S)
            S = S(:);
            delta_S = obj.model_fun(S);
            if obj.t == 1; Sold = obj.S0(:); else; Sold = obj.stores(obj.t-1,:)'; end
            err = (S - Sold)/obj.delta_t - delta_S';
        end 
        
        % SOLVE_STORES solves the stores ODEs 
        function [Snew, resnorm, solver, iter] = solve_stores(obj, Sold)
            
            solver_opts = obj.solver_opts;
            
            % This reduces the tolerance to a fraction of the smallest store,
            % if stores are very small, with 1E-6 as minimum
            % (if resnorm_tolerance is 0.1 as default)
            resnorm_tolerance = solver_opts.resnorm_tolerance * min(min(Sold) + 1E-5, 1);
            
            % create vectors for each of the three solutions (NewtonRaphon,
            % fsolve and lsqnonlin), this way if all three run it takes the
            % best at the end and not the last one.
            Snew_v    = zeros(3, obj.numStores);
            resnorm_v = Inf(3, 1);
            iter_v    = ones(3,1);
            
            % first try to solve the ODEs using NewtonRaphson
            [tmp_Snew, tmp_fval] = ...
                            NewtonRaphson(@obj.solve_fun_IE,...
                                          Sold,...
                                          solver_opts.NewtonRaphson);
            tmp_resnorm = sum(tmp_fval.^2);
            
            Snew_v(1,:)  = tmp_Snew;
            resnorm_v(1) = tmp_resnorm;
            
            % if NewtonRaphson doesn't find a good enough solution, run FSOLVE
            if tmp_resnorm > resnorm_tolerance  
                [tmp_Snew,tmp_fval,~,tmp_iter] = ...
                            rerunSolver('fsolve', ...              
                                        solver_opts.fsolve, ...            % solver options
                                        @obj.solve_fun_IE,...              % system of ODEs
                                        solver_opts.rerun_maxiter, ...     % maximum number of re-runs
                                        resnorm_tolerance, ...             % convergence tolerance
                                        tmp_Snew, ...                      % recent estimates
                                        Sold, ...                          % storages at previous time step
                                        obj.store_min, ...                 % lower bounds
                                        obj.store_max);                    % upper bounds 
                
                tmp_resnorm = sum(tmp_fval.^2);
                
                Snew_v(2,:)  = tmp_Snew;
                resnorm_v(2) = tmp_resnorm;
                iter_v(2)    = tmp_iter;

                % if FSOLVE doesn't find a good enough solution, run LSQNONLIN
                if tmp_resnorm > resnorm_tolerance
                    [tmp_Snew,tmp_fval,~,tmp_iter] = ...
                            rerunSolver('lsqnonlin', ...              
                                        solver_opts.lsqnonlin, ...         % solver options
                                        @obj.solve_fun_IE,...              % system of ODEs
                                        solver_opts.rerun_maxiter, ...     % maximum number of re-runs
                                        resnorm_tolerance, ...             % convergence tolerance
                                        tmp_Snew, ...                      % recent estimates
                                        Sold, ...                          % storages at previous time step
                                        obj.store_min, ...                 % lower bounds
                                        obj.store_max);                    % upper bounds
                    
                    tmp_resnorm = sum(tmp_fval.^2);

                    Snew_v(3,:)  = tmp_Snew;
                    resnorm_v(3) = tmp_resnorm;
                    iter_v(3)    = tmp_iter;
                    
                end
            end
            
            % get the best solution
            [resnorm, solver_id] = min(resnorm_v);
            Snew = Snew_v(solver_id,:);
            iter = iter_v(solver_id);
            solvers = ["NewtonRaphson", "fsolve", "lsqnonlin"];
            solver = solvers(solver_id);
            
        end
        
        % RUN runs the model with a given climate input, initial stores,
        % parameter set and solver settings.
        % none of the arguments are needed, they can be set beforehand with
        % obj.theta = theta; obj.input_climate = input_climate; etc.
        % then simply obj.run() without arguments.
        function [fluxes, stores] = run(obj,...
                                        input_climate,...         
                                        S0,...
                                        theta,...
                                        solver_opts)

            if nargin > 4 && ~isempty(solver_opts)
                obj.solver_opts = solver_opts;
            end
            if nargin > 3 && ~isempty(theta)
                obj.theta = theta;
            end
            if nargin > 2 && ~isempty(S0)
                obj.S0 = S0;
            end
            if nargin > 1 && ~isempty(input_climate)
                obj.input_climate = input_climate;
            end
            
            
            % run INIT_ method, this will calculate all auxiliary parameters
            % and set up routing vectors and store limits
            obj.init_();
            
            % if the max span of the stores is small (less than 2 orders of
            % magnetude as the tolerance), reduce the tolerance
%             if min(obj.store_max - obj.store_min) < obj.solver_opts.resnorm_tolerance*100
%                 obj.solver_opts.resnorm_tolerance = obj.solver_opts.resnorm_tolerance/100;
%             end

            t_end = size(obj.input_climate, 1);
            
            for t = 1:t_end
               obj.t = t;
               if t == 1; Sold = obj.S0(:); else; Sold = obj.stores(t-1,:)'; end
               %obj.input_climate = [P(t) Ep(t) T(t)];
               [Snew,resnorm,solver,iter] = obj.solve_stores(Sold);
               
               [dS, f] = obj.model_fun(Snew);
    
               obj.fluxes(t,:) = f * obj.delta_t;
               obj.stores(t,:) = Sold + dS' * obj.delta_t;
               
               obj.step();
               
               obj.solver_data.resnorm(t) = resnorm;
               obj.solver_data.solver(t) = solver;
               obj.solver_data.iter(t) = iter;
               
            end
            fluxes = obj.fluxes;
            stores = obj.stores;
            obj.status = 1;
            
        end
        
        % GET_OUTPUT runs the model exactly like RUN, but output is
        % consistent with current MARRMoT
        function [fluxOutput,...
                  fluxInternal,...
                  storeInternal,...
                  waterBalance,...
                  solverSteps] = get_output(obj,...
                                            varargin)
            
            if nargin > 1 || isempty(obj.status) || obj.status == 0 
                [~] = obj.run(varargin{:});
            end
            
            % --- Fluxes leaving the model ---
            fluxOutput.Ea     = sum(obj.fluxes(:,obj.FluxGroups.Ea),2)';
            fluxOutput.Q      = sum(obj.fluxes(:,obj.FluxGroups.Q),2)';
            
            % --- Fluxes internal to the model ---
            fluxInternal = struct;
            for i = 1:obj.numFluxes
                fluxInternal.(obj.FluxNames(i)) = obj.fluxes(:,i)';
            end
            
            % --- Stores ---
            storeInternal = struct;
            for i = 1:obj.numStores
                storeInternal.(obj.StoreNames(i)) = obj.stores(:,i)';
            end
            
            % --- Water balance, if requested ---
            if nargout >= 4
                waterBalance = obj.check_waterbalance();
            end
            
            % --- step-by-step data of the solver, if requested ---
            if nargout == 5
                solverSteps = obj.solver_data;
            end
        end
        
        % CHECK_WATERBALANCE returns the waterbalance given fluxes and
        % stores (i.e. the output of RUN)
        % like in MARRMoT1, it will print to screen
        function [out] = check_waterbalance(obj)
            
            if isempty(obj.status) || obj.status == 0 
                error('run model with obj.run() before requesting the water balance');
            end
            
            % Get variables
            P  = obj.input_climate(:,1);
            fg = fieldnames(obj.FluxGroups);
            OutFluxes = zeros(1,numel(fg));
            for k=1:numel(fg)                                              % cumulative of each flow leaving the model
                idx = abs(obj.FluxGroups.(fg{k}));
                signs = sign(obj.FluxGroups.(fg{k}));
                OutFluxes(k) = sum(sum(signs.*obj.fluxes(:,idx), 1),2);
            end
            if isempty(obj.StoreSigns); obj.StoreSigns = repelem(1, obj.numStores); end
            dS = obj.StoreSigns .* (obj.stores(end,:) - obj.S0');          % difference of final and initial storage for each store
            if isempty(obj.fluxes_stf); obj.fluxes_stf = {}; end
            R = cellfun(@sum, obj.fluxes_stf);                             % cumulative of each flows still to be routed
            
            % calculate water balance
            out = sum(P) - ...                                             % input from precipitation
                sum(OutFluxes) - ...                                       % all fluxes leaving the model (some may be entering, but they should have a negative sign)
                sum(dS) - ...                                              % all differences in storage
                sum(R);                                                    % all flows still being routed
            
            disp(['Total P  = ',num2str(sum(P)),' mm.'])
            for k = 1:numel(fg)
                disp(['Total ',char(fg(k)),' = ',...
                      num2str(-OutFluxes(k)),' mm.'])
            end
            for s = 1:obj.numStores
                if obj.StoreSigns(s) == -1
                    ending=' (deficit store).';
                else
                    ending='.';
                end
                disp(['Delta S',num2str(s),' = ',...
                      num2str(-dS(s)),' mm',ending])
            end
            if ~isempty(R)
                disp(['On route = ',num2str(-sum(R)),' mm.'])
            end

        disp('-------------')
        disp(['Water balance = ', num2str(out), ' mm.'])
        end
        
        % GET_STREAMFLOW only returns the streamflow, runs the model if it
        % hadn't run already.
        function Q = get_streamflow(obj,...
                                    varargin)
            
            if nargin > 1 || isempty(obj.status) || obj.status == 0
                [~] = obj.run(varargin{:});
            end
        
            Q = sum(obj.fluxes(:,obj.FluxGroups.Q),2);
        end
        
        % CALC_PAR_FITNESS calculates the fitness of the set of parameters
        % in obj.theta, needs the model to have run already
        function fitness = calc_fitness(obj,...
                                        Q_obs,...                          % observed streamflow
                                        of_name,...                        % name of the objective function
                                        fit_idx,...                        % indices to use to calculate the fitness function
                                        varargin)                          % additional arguments to objective function
            
            if isempty(obj.status) || obj.status == 0 
                error('run model with obj.run() before requesting fitness');
            end
                                        
            Q_sim = sum(obj.fluxes(:,obj.FluxGroups.Q),2);
            fitness = feval(of_name, Q_obs, Q_sim, fit_idx, varargin{:});
        end
        
        % CALIBRATE uses the chosen algorithm to find the optimal parameter
        % set, given model inputs, objective function and observed streamflow.
        % the function chosen in algorithm should have the same inputs and
        % outputs as MATLAB's fminsearch.
        function  [par_opt,...                                             % optimal parameter set
                   of_cal,...                                              % value of objective function at par_opt
                   stopflag,...                                            % flag indicating reason the algorithm stopped
                   output] = ...                                           % output, see fminsearch for detail
                             calibrate(obj,...
                                       Q_obs,...                           % observed streamflow
                                       cal_idx,...                         % timesteps to use for model calibration
                                       optim_fun,...                       % function to use for optimisation (must have same structure as fminsearch)
                                       par_ini,...                         % initial parameter estimates
                                       optim_opts,...                      % options to optim_fun
                                       of_name,...                         % name of objective function to use
                                       inverse_flag,...                    % should the OF be inversed?
                                       varargin)                           % additional arguments to the objective function
             
             if isempty(obj.input_climate) || isempty(obj.delta_t) ||...
                     isempty(obj.S0) || isempty(obj.solver_opts)
                 error(['input_climate, delta_t, S0 and solver_opts '...
                        'attributes must be specified before calling '...
                        'calibrate. Use obj.input_climate = ... etc.']);
             end
             
             % if the list of timesteps to use for calibration is empty,
             % use all steps 
             if isempty(cal_idx)
                 cal_idx = 1:length(Q_obs);
             end
                                   
             % helper function to calculate fitness given a set of
             % parameters
             function fitness = fitness_fun(par)
                 [~] = obj.run([],[],par);
                 Q_sim = sum(obj.fluxes(:,obj.FluxGroups.Q),2);
                 fitness = (-1)^inverse_flag*feval(of_name, Q_obs, Q_sim, cal_idx, varargin{:});
             end
             
             % if the initial parameter set isn't set,  start from mean
             % values of parameter range
             if isempty(par_ini)
                 %par_ini = obj.parRanges(:,1) + rand(obj.numParams,1).*(obj.parRanges(:,2)-obj.parRanges(:,1));
                 par_ini = mean(obj.parRanges,2);
             end
             
             [par_opt,...                                                  % optimal parameter set at the end of the optimisation
                 of_cal,...                                                % value of the objective function at par_opt
                 stopflag,...                                              % flag indicating reason the algorithm stopped
                 output] = ...                                             % output, see fminsearch for detail
                           feval(optim_fun,...                             % run the optimisation algorithm chosen
                                 @fitness_fun,...                          % function to optimise is the fitness function
                                 par_ini,...                               % initial parameter set
                                 optim_opts);                              % optimiser options
             
             % if of_cal was inverted, invert it back before returning
             of_cal = (-1)^inverse_flag * of_cal;
        end
         
         % function to return default solver options
         function solver_opts = default_solver_opts(obj)
            solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
            solver_opts.rerun_maxiter   = 6;                                           % Maximum number of re-runs
            solver_opts.NewtonRaphson = optimset('MaxIter', obj.numStores * 10);
            solver_opts.fsolve = optimoptions('fsolve',...
                                              'Display','none',...                     % Disable display settings
                                              'JacobPattern', obj.JacobPattern);
            solver_opts.lsqnonlin = optimoptions('lsqnonlin',...                       % lsqnonlin settings for cases where fsolve fails
                                                 'Display','none',...
                                                 'JacobPattern',obj.JacobPattern,...
                                                 'MaxFunEvals',1000);
         end
         
         % function to add new solver opts to the default ones
         function solver_opts = add_to_def_opts(obj, opts)
             def_opts = obj.default_solver_opts();             
             if nargin == 1 || isempty(opts)
                 solver_opts = def_opts;
             else
                 def_fields = fields(def_opts);
                 % for each field in the default options (5 at the moment)
                 for k = 1:length(def_fields)
                     field = def_fields{k};
                     % if the field is not provided, use the default one
                     if ~isfield(opts, field) || isempty(opts.(field))
                         solver_opts.(field) = def_opts.(field);
                     % if the field is provided, and the dafault is a struct,
                     % add the new values to the default struct
                     elseif isstruct(def_opts.(field))
                         solver_opts.(field) = optimset(def_opts.(field),opts.(field));
                     % if the field is provided, and the dafault is not a struct,
                     % discard the default and use the new value
                     else
                         solver_opts.(field) = opts.(field);
                     end
                 end
             end             
         end
    end
end
        
