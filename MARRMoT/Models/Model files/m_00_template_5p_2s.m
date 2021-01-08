classdef m_00_template_5p_2s < MARRMoT_model
    % Class for MARRMoT model teplate
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_00_template_5p_2s(delta_t, theta)
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 6;                                             % number of model fluxes
            obj.numParams = 5;
            
            obj.JacobPattern  = [1 1;...
                                 1 1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 40;      % Smax [mm]
                             0   , 2 ;      % kc, capillary rise [mm/d]
                             0   , 3 ;      % kp, percolation rate [mm/d]
                             0.5 , 1 ;      % ks, base flow time parameter [d-1]
                             1   , 5];      % time delay of routing scheme [d]
            
            obj.StoreNames = ["S1" "S2"];                                  % Names for the stores
            obj.FluxNames  = ["cap", "ea", "qo", "perc", "qs","qt"];       % Names for the fluxes
            
            obj.Flux_Ea_idx = 2;                                           % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = 6;                                           % Index or indices of fluxes to add to Streamflow
            
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
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            delay   = theta(5);     % Routing delay [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(delay,delta_t);
            
            obj.uhs        = {uh};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            S1max   = theta(1);     % Maximum soil moisture storage [mm]
            kc      = theta(2);     % Maximum capillary rise [mm/d]
            kp      = theta(3);     % Maximum percolation [mm/d]
            ks      = theta(4);     % Runoff coefficient [d-1]
            delay   = theta(5);     % Routing delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            
            % climate input
            c = obj.input_climate;
            P  = c(1);
            Ep = c(2);
            T  = c(3);
            
            % fluxes functions
            flux_cap  = min(max(kc*(S1max-S1)/S1max,0),S2/delta_t);
            flux_ea   = min(S1/delta_t,Ep);
            flux_qo   = P.*(1-smoothThreshold_storage_logistic(S1,S1max,0.001));
            flux_perc = min(kp.*S1/S1max,S1/delta_t);
            flux_qs   = ks*S2;
            flux_qt   = uh(1) * (flux_qo + flux_qs) + stf(1);
            

            % stores ODEs
            dS1 = P         + flux_cap - flux_ea  - flux_qo - flux_perc;
            dS2 = flux_perc - flux_qs  - flux_cap;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_cap,  flux_ea, flux_qo,...
                      flux_perc, flux_qs, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % input fluxes to the unit hydrographs 
            flux_qo = fluxes(3);
            flux_qs = fluxes(5);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf      = (uh .* (flux_qo + flux_qs)) + stf;
            stf      = circshift(stf,-1);
            stf(end) = 0;
            
            obj.fluxes_stf = {stf};
        end
    end
end