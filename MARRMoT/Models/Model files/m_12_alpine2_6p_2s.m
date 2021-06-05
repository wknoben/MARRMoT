classdef m_12_alpine2_6p_2s < MARRMoT_model
    % Class for alpine2 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_12_alpine2_6p_2s()
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
            
        end
        
        % INITialisation function
        function obj = init(obj)
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
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
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
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end