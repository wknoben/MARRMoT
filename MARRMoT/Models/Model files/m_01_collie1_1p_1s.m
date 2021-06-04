classdef m_01_collie1_1p_1s < MARRMoT_model
    % Class for collie1 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_01_collie1_1p_1s(delta_t, theta)
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 2;                                             % number of model fluxes
            obj.numParams = 1;
            
            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [1   , 2000];      % Smax [mm]
            
            obj.StoreNames = ["S1"];                                       % Names for the stores
            obj.FluxNames  = ["ea", "qse"];                                % Names for the fluxes
            
            obj.FluxGroups.Ea = 1;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 2;                                         % Index or indices of fluxes to add to Streamflow
            
            % setting delta_t and theta triggers the function obj.init()
            if nargin > 0 && ~isempty(delta_t)
                obj.delta_t = delta_t;
            end
            if nargin > 1 && ~isempty(theta)
                obj.theta = theta;
            end
        end
        
        % INIT is run automatically as soon as both theta and delta_t are
        % set (it is therefore ran only once at the beginning of the run. 
        % Use it to initialise all the model parameters (in case there are
        % derived parameters) and unit hydrographs and set minima and
        % maxima for stores based on parameters.
        function obj = init(obj)
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            S1max   = theta(1);                         % Maximum soil moisture storage [mm]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);

            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            
            % fluxes functions
            flux_ea   = evap_7(S1, S1max, Ep, delta_t);
            flux_qse  = saturation_1(P,S1,S1max);
            
            % stores ODEs
            dS1 = P - flux_ea  - flux_qse;
            
            % outputs
            dS = [dS1];
            fluxes = [flux_ea, flux_qse];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end