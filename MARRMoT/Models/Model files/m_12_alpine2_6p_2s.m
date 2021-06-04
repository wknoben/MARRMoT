classdef m_12_alpine2_6p_2s < MARRMoT_model
    % Class for alpine2 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_12_alpine2_6p_2s(delta_t, theta)
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 7;                                             % number of model fluxes
            obj.numParams = 6; 

            obj.JacobPattern  = [1,0;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [-3  , 5        % tt [celsius]
                             0   , 20;      % degree-day-factor [mm/degree celsius/d]
                             1   , 2000;    % Smax [mm]
                             0.05, 0.95;    % Field capacity as fraction of Smax [-]
                             0   , 1;       % time delay of interflow [d-1]
                             0   , 1];      % time delay of baseflow [d-1]
            
            obj.StoreNames = ["S1" "S2"];                                  % Names for the stores
            obj.FluxNames  = ["ps", "pr", "qn", "ea", "qse", "qin", "qbf"];% Names for the fluxes
            
            obj.FluxGroups.Ea = [4];                                       % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [5 6 7];                                   % Index or indices of fluxes to add to Streamflow
            
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
            tt   = theta(1);     % Threshold temperature for snowfall/snowmelt [celsius]
            ddf  = theta(2);     % Degree-day-factor for snowmelt [mm/d/celsius]
            Smax = theta(3);     % Maximum soil moisture storage [mm]
            Cfc  = theta(4);     % Field capacity as fraction of Smax [-]
            tcin = theta(5);     % Interflow runoff coefficient [d-1]
            tcbf = theta(6);     % Baseflow runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ps  = snowfall_1(P,T,tt);
            flux_pr  = rainfall_1(P,T,tt);
            flux_qn  = melt_1(ddf,tt,T,S1,delta_t);
            flux_ea  = evap_1(S2,Ep,delta_t);
            flux_qse = saturation_1(P+flux_qn,S2,Smax);
            flux_qin = interflow_8(S2,tcin,Cfc*Smax);
            flux_qbf = baseflow_1(tcbf,S2);

            % stores ODEs
            dS1 = flux_ps - flux_qn;
            dS2 = flux_pr + flux_qn - flux_ea - flux_qse - flux_qin - flux_qbf;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_ps, flux_pr, flux_qn, flux_ea,...
                      flux_qse, flux_qin, flux_qbf];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end