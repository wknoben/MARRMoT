classdef m_nn_template_pp_ss < MARRMoT_model
% Class for hydrologic conceptual model: template

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_nn_template_pp_ss()
            obj.numStores = s;                                             % number of model stores
            obj.numFluxes = f;                                             % number of model fluxes
            obj.numParams = p;
            
            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
                         
            %                min, max
            obj.parRanges = [0  , Inf;      % parameter 1 [unit]
                             0  , Inf;      % parameter 2 [unit]
                             
                             0  , Inf];     % parameter p [unit]
            
            obj.StoreNames = {"S1", "S2", "Ss"};                            % Names for the stores
            obj.FluxNames  = {"f1", "f2", "f2", "r1", "r2"};               % Names for the fluxes
            
            obj.FluxGroups.Ea = 0;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 0;                                         % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INIT is called by INIT_ and runs once beofore each model run. 
        % Use it to initialise all the model parameters (in case there are
        % derived parameters) and unit hydrographs and set minima and
        % maxima for stores based on parameters.
        function obj = init(obj)
            % extract theta and delta_t from attributes
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            % needed parameters
            p1 = theta(1); % parameter 1 [unit]
            p2 = theta(2); % parameter 2 [unit]
            p5 = theta(5); % parameter 5 [unit]
            p6 = theta(6); % parameter 6 [unit]
            
            % min and max of stores
            obj.store_mim(1) = p1;                      % maxima and minima of individual stores
            obj.store_max    = 10 .* [p1, p2, p1+p2];   % or as a whole
            
            % unit hydrographs         
            uh1 = uh1_function(p5,delta_t);
            uh2 = uh2_function(p6,delta_t);
            obj.uhs = {uh1 uh2}; 
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            p1 = theta(1); % parameter 1 [unit]
            p2 = theta(2); % parameter 2 [unit]
            p3 = theta(3); % parameter 3 [unit]
            p4 = theta(4); % parameter 4 [unit]
            pp = theta(p); % parameter p [unit]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs
            uhs = obj.uhs;
            uh1 = uhs{1};
            uh2 = uhs{2};
            
            % stores (S is the input to the method, not an attribute)
            S1 = S(1);
            S2 = S(2);
            Ss = S(s);
            
            % climate input
            t = obj.t;                    % this time step
            c = obj.input_climate(t,:);   % climate at this step
            P  = c(1);
            Ep = c(2);
            T  = c(3);
            
            % fluxes functions
              % use flux functions from flux files
            flux_f1  = flux_function1(arguments);   % arguments of flux functions are 
            flux_f2  = flux_function2(arguments);   % parameters, store values, climate inputs
            flux_f3  = flux_function3(arguments);   % or other fluxes
              % or the function route(flux_in, uh) to rout a flux through a
              % unit hydrograph
            flux_r1 = route(flux_f3, uh1);
            flux_r2 = route(0.5 .* (flux_f1 + flux_r1), uh2);
            
            % stores ODEs
            dS1 = P       + flux_f1 - flux_r1;         % each store has 1 ODE,                                    
            dS2 = flux_f2 - flux_f3 - flux_r2;         % in which entering fluxes are added
            dSs = enternig_fluxes - exiting_fluxes;    % and exiting_fluxes are subtracted
            
            % outputs
            dS = [dS1 dS2 dSs];                       % output are arrays of all stores dS and
            fluxes = [flux_f1,  flux_f2, flux_f3,...  % all fluxes at this timestep. Order must match
                      flux_r1, flux_r2};              % the naming in obj.StoreNames and obj.FluxNames
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh1 = uhs{1};
            uh2 = uhs{2};
            
            % input fluxes to the unit hydrographs at this timestep 
            fluxes = obj.fluxes(obj.t,:); 
            flux_f1 = fluxes(1);
            flux_f3 = fluxes(3);
            flux_r1 = fluxes(4);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh1 = update_uh(uh1, flux_f3);
            uh2 = update_uh(uh2, 0.5 .* (flux_f1 + flux_r1));
            
            obj.uhs = {uh1, uh2};
        end
    end
end