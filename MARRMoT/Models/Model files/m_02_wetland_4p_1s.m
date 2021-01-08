classdef m_02_wetland_4p_1s < MARRMoT_model
    % Class for wetland model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_02_wetland_4p_1s(delta_t, theta)
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 5;                                             % number of model fluxes
            obj.numParams = 4;

            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [0   , 5;      % Dw, interception capacity [mm]
                             0   , 10;     % Betaw, soil misture distribution parameter [-]
                             1   , 2000;   % Swmax, soil misture depth [mm]
                             0   , 1];     % kw, base flow time parameter [d-1]
            
            obj.StoreNames = ["S1"];                                       % Names for the stores
            obj.FluxNames  = ["pe",  "ei",  "ew",  "qwsof", "qwgw"];       % Names for the fluxes
            
            obj.Flux_Ea_idx = [2 3];                                       % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = [4 5];                                       % Index or indices of fluxes to add to Streamflow
            
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
            dw    = theta(1);     % Daily interception [mm]
            betaw = theta(2);     % Soil moisture storage distribution parameter [-]
            swmax = theta(3);     % Maximum soil moisture storage [mm], 
            kw    = theta(4);     % Runoff coefficient [d-1]
            
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pe    = interception_2(P,dw);
            flux_ei    = P - flux_pe;
            flux_ew    = evap_1(S1,Ep,delta_t);
            flux_qwsof = saturation_2(S1,swmax,betaw,flux_pe);
            flux_qwgw  = baseflow_1(kw,S1);

            % stores ODEs
            dS1 = flux_pe - flux_ew - flux_qwsof - flux_qwgw; 
            
            % outputs
            dS = [dS1];
            fluxes = [flux_pe,  flux_ei,  flux_ew,  flux_qwsof, flux_qwgw];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end