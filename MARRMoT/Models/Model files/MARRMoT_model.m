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
            err = (S - obj.Sold')/obj.delta_t - delta_S';
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
            
            for t = 1:t_end
               obj.t = t;
               if t == 1; obj.Sold = S0; else; obj.Sold = stores(t-1,:); end
               obj.input_climate = [P(t) Ep(t) T(t)];
               Snew = obj.solve_stores(solver_opts);
               
               [dS, f] = obj.model_fun(Snew);
    
               fluxes(t,:) = f * obj.delta_t;
               stores(t,:) = obj.Sold + dS * obj.delta_t;
               
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
                                              solver)
            
            [fluxes, stores] = obj.run(fluxInput, storeInitial, solver);
            
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
    end
end
        
