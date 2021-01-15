classdef MARRMoT_model < handle
    % Class for MARRMoT models
    properties
        numStores         % number of model stores
        numFluxes         % number of model fluxes
        numParams         % number of model parameters
        parRanges         % default parameter ranges
        JacobPattern      % pattern of the Jacobian matrix of model store ODEs
        StoreNames        % Names for the stores
        FluxNames         % Names for the fluxes
        Flux_Ea_idx       % indices of fluxes related to Actual ET
        Flux_Q_idx        % indices of fluxes related to Streamflow
        theta             % Set of parameters
        delta_t           % time step
        store_min         % store minimum values
        store_max         % store maximum values
        uhs               % unit hydrographs
        fluxes_stf        % still-to-flow fluxes
        t                 % current timestep
        Sold              % Vector of stores at time t-1
        input_climate     % Vector of climate at t
    end
    methods
        
        % Set methods for delta_t and theta are set up so that when both
        % delta_t and theta are set, the method INIT is called, allowing to
        % initiate model characteristics based on theta and delta_t (namely
        % unit hydrographs, auxiliary parameters or store limits).
        function [] = set.delta_t(obj, value)
            if length(value) == 1
                obj.delta_t = value;
                if ~isempty(obj.theta); obj.init(); end
            else
                error('delta_t must be a scalar')
            end
        end
        function [] = set.theta(obj, value)
            if length(value) == obj.numParams
                obj.theta = value;
                if ~isempty(obj.delta_t); obj.init(); end
            else
                error(['theta must have ' int2str(obj.numParams) ' elements'])
            end
        end
        
        % ODE approximation with Implicit Euler time-stepping scheme
        function err = solve_fun_IE(obj, S)
            S = S(:);
            delta_S = obj.model_fun(S);
            err = (S - obj.Sold)/obj.delta_t - delta_S';
        end 
        
        % SOLVE_STORES solves the stores ODEs 
        function Snew = solve_stores(obj, solver_opts)

            [Snew, fval, exitflag] = NewtonRaphson(@obj.solve_fun_IE,...
                                                   obj.Sold,...
                                                   solver_opts.NewtonRaphson);
            resnorm = sum(fval.^2);
            
            % if NewtonRaphson doesn't find a good enough solution, run FSOLVE
            if exitflag <= 0 || resnorm > solver_opts.resnorm_tolerance
                [Snew,fval,exitflag] = fsolve(@obj.solve_fun_IE,...        % system of storage equations
                                              Snew,...                     % storage values at the end of NewtonRaphson solver
                                              solver_opts.fsolve);         % solver options
                resnorm = sum(fval.^2);
                
                % if FSOLVE doesn't find a good enough solution, run LSQNONLIN
                if exitflag <= 0 || resnorm > solver_opts.resnorm_tolerance
                    [Snew,~,~] = rerunSolver('lsqnonlin', ...              
                        solver_opts.lsqnonlin, ...                         % solver options
                        @obj.solve_fun_IE,...                              % system of ODEs
                        solver_opts.resnorm_maxiter, ...                   % maximum number of re-runs
                        solver_opts.resnorm_tolerance, ...                 % convergence tolerance
                        Snew, ...                                          % recent estimates
                        obj.Sold, ...                                      % storages at previous time step
                        obj.store_min, ...                                 % lower bounds
                        obj.store_max);                                    % upper bounds 
                end
            end
        end
        
        % RUN runs the model with a given climate input, initial stores,
        % parameter set and solver settings.
        function [fluxes, stores] = run(obj,...
                                        fluxInput,...
                                        storeInitial,...
                                        solver_opts)

            P     = fluxInput.precip./obj.delta_t;          % [mm/delta_t] -> [mm/d]       
            Ep    = fluxInput.pet./obj.delta_t;             % [mm/delta_t] -> [mm/d]       
            T     = fluxInput.temp;
            t_end = length(P);
            
            S0 = storeInitial;
            
            stores = zeros(t_end, obj.numStores);
            fluxes = zeros(t_end, obj.numFluxes);
            
            % use default options for the solvers if they are not set
            if nargin < 4 || isempty(solver_opts)
                solver_opts = obj.default_solver_opts;
            end
            
            for t = 1:t_end
               obj.t = t;
               if t == 1; obj.Sold = S0(:); else; obj.Sold = stores(t-1,:)'; end
               obj.input_climate = [P(t) Ep(t) T(t)];
               Snew = obj.solve_stores(solver_opts);
               
               [dS, f] = obj.model_fun(Snew);
    
               fluxes(t,:) = f * obj.delta_t;
               stores(t,:) = obj.Sold + dS' * obj.delta_t;
               
               obj.step(f);
            end
        end
        % GET_OUTPUT runs the model exactly like RUN, but output is
        % consistent with current MARRMoT
        function [fluxOutput,...
                  fluxInternal,...
                  storeInternal] = get_output(obj,...
                                              fluxInput,...
                                              storeInitial,...
                                              solver_opts)
            
            [fluxes, stores] = obj.run(fluxInput, storeInitial, solver_opts);
            
            % --- Fluxes leaving the model ---
            fluxOutput.Ea     = sum(fluxes(:,obj.Flux_Ea_idx),2)';
            fluxOutput.Q      = sum(fluxes(:,obj.Flux_Q_idx),2)';
            
            % --- Fluxes internal to the model ---
            fluxInternal = struct;
            for i = 1:obj.numFluxes
                fluxInternal.(obj.FluxNames(i)) = fluxes(:,i)';
            end
            
            % --- Stores ---
            storeInternal = struct;
            for i = 1:obj.numStores
                storeInternal.(obj.StoreNames(i)) = stores(:,i)';
            end
        end
        
        % GET_STREAMFLOW runs the model and only returns the streamflow, it
        % is useful for calibration
        function Q = get_streamflow(obj,...
                                    fluxInput,...                  % struct of input fluxes
                                    storeInitial,...               % initial values of stores
                                    solver_opts)
             fluxes = obj.run(fluxInput, storeInitial, solver_opts);
             Q = sum(fluxes(:,obj.Flux_Q_idx),2);
        end
        
        % CALC_PAR_FITNESS calculates the fitness of the set of parameters
        % in obj.theta based on a given model inputs, objective function
        % and observed streamflow
        function fitness = calc_par_fitness(obj,...
                                            fluxInput,...                  % struct of input fluxes
                                            storeInitial,...               % initial values of stores
                                            solver_opts,...                % solver options
                                            Q_obs,...                      % observed streamflow
                                            of_name,...                    % name of the objective function
                                            varargin)                      % additional arguments to
            Q_sim = obj.get_streamflow(fluxInput, storeInitial, solver_opts);
            fitness = feval(of_name, Q_obs, Q_sim, varargin{:});
        end
        
        % CALIBRATE uses CMA-ES to find the optimal parameter set given
        % model inputs, objective function and observed streamflow                                           
         function [par_opt,...                                             % optimal parameter set
                   of_cal,...                                              % value of objective function at calibration
                   n_feval,...                                             % number of function evaluations
                   stopflag]= calibrate(obj,...
                                        fluxInput,...                      % struct of input fluxes
                                        storeInitial,...                   % initial values of stores
                                        solver_opts,...                    % solver options
                                        Q_obs,...                          % observed streamflow
                                        cmaes_sigma0,...                   % initial sigma for CMA-ES
                                        cmaes_opts,...                     % CMA-ES options
                                        of_name,...                        % name of objective function to use
                                        inverse_flag,...                   % should the OF be inversed?
                                        varargin)                          % additional arguments to
             
             % helper function to calculate fitness given a set of
             % parameters
             function fitness = fitness_fun(par)
                 obj.theta = par;
                 Q_sim = obj.get_streamflow(fluxInput, storeInitial, solver_opts);
                 fitness = (-1)^inverse_flag*feval(of_name, Q_obs, Q_sim, varargin{:});
             end
             
             % start from mean parameter set
             %par_ini = obj.parRanges(:,1) + rand(obj.numParams,1).*(obj.parRanges(:,2)-obj.parRanges(:,1));
             par_ini = mean(obj.parRanges,2);
             
             % set default options for cmaes if options are empty
             [def_cmaes_opts, def_cmaes_sigma0] = obj.default_cmaes_opts();
             if isempty(cmaes_opts);   cmaes_opts   = def_cmaes_opts;   end
             if isempty(cmaes_sigma0); cmaes_sigma0 = def_cmaes_sigma0; end
             
             % set parameter bounds if they are not set
             if ~isfield(cmaes_opts, 'LBounds') || isempty(cmaes_opts.LBounds)
                 cmaes_opts.LBounds = obj.parRanges(:,1);
             end
             if ~isfield(cmaes_opts, 'UBounds') || isempty(cmaes_opts.UBounds)
                 cmaes_opts.UBounds = obj.parRanges(:,2);
             end
             
             % run CMA-ES
             [par_opt,...
                 of_cal,...
                 n_feval,...
                 stopflag] = cmaes(@fitness_fun,...
                                   par_ini,...
                                   cmaes_sigma0,...
                                   cmaes_opts);
         end
         
         % function to return default solver options
         function solver_opts = default_solver_opts(obj)
            solver_opts.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
            solver_opts.resnorm_maxiter   = 6;                                         % Maximum number of re-runs
            solver_opts.NewtonRaphson = optimset('MaxIter', obj.numStores * 10);
            solver_opts.fsolve = optimoptions('fsolve',...
                                              'Display','none',...                     % Disable display settings
                                              'JacobPattern', obj.JacobPattern);
            solver_opts.lsqnonlin = optimoptions('lsqnonlin',...                       % lsqnonlin settings for cases where fsolve fails
                                                 'Display','none',...
                                                 'JacobPattern',obj.JacobPattern,...
                                                 'MaxFunEvals',1000);
         end
         
         % function to return default CMA-ES options
         function [cmaes_opts, cmaes_sigma0] = default_cmaes_opts(obj)
            cmaes_opts.LBounds  = obj.parRanges(:,1);                      % lower bounds of parameters
            cmaes_opts.UBounds  = obj.parRanges(:,2);                      % upper bounds of parameters
            cmaes_opts.PopSize  = 2 * (4 + floor(3*log(obj.numParams)));   % population size is 2x the defaul
            
            % starting sigma
            cmaes_sigma0 = .3*(obj.parRanges(:,2) - obj.parRanges(:,1));   % starting sigma (this is default)
            
            % stopping criteria
            cmaes_opts.TolX       = 1e-6 * min(cmaes_sigma0);              % stopping criterion on changes to parameters 
            cmaes_opts.TolFun     = 1e-4;                                  % stopping criterion on changes to fitness function
            cmaes_opts.TolHistFun = 1e-5;                                  % stopping criterion on changes to fitness function
         end
    end
end
        
