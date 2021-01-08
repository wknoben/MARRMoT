classdef m_14_topmodel_7p_2s < MARRMoT_model
    % Class for topmodel model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_14_topmodel_7p_2s(delta_t, theta)
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 6;                                             % number of model fluxes
            obj.numParams = 7;
            
            obj.JacobPattern  = [1 1;...
                                 1 1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1, 2000;    % suzmax, Maximum soil moisture storage in unsatured zone [mm]
                             0.05, 0.95; % st, Threshold for flow generation and evap change as fraction of suzmax [-]
                             0,1;        % kd, Leakage to saturated zone flow coefficient [mm/d]
                             0.1,200;    % q0, Zero deficit base flow speed [mm/d]
                             0,1;        % m, Baseflow coefficient [mm-1]
                             1, 7.5;     % chi, Gamma distribution parameter [-]
                             0.1, 5];    % phi, Gamma distribution parameter [-]
            
            obj.StoreNames = ["S1" "S2"];                                  % Names for the stores
            obj.FluxNames  = ["qof", "peff", "ea", "qex", "qv", "qb"];     % Names for the fluxes
            
            obj.Flux_Ea_idx = 3;                                           % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = [1 4 6];                                     % Index or indices of fluxes to add to Streamflow
            
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
            suzmax  = theta(1);     % Maximum soil moisture storage in unsatured zone [mm]
            st      = theta(2);     % Threshold for flow generation and evap change as fraction of suzmax [-]
            kd      = theta(3);     % Leakage to saturated zone flow coefficient [mm/d]
            q0      = theta(4);     % Zero deficit base flow speed [mm/d]
            f       = theta(5);     % Baseflow scaling coefficient [mm-1]
            chi     = theta(6);     % Gamma distribution parameter [-]
            phi     = theta(7);     % Gamma distribution parameter [-]
            mu      = 3;            % Gamma distribution parameter, fixed (Clark et al, 2008)
            lambda  = chi*phi+mu;   % Ac computation parameter, mean of the gamma distribution

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
            flux_qof  = saturation_7(chi,phi,3,lambda,f,S2,P);
            flux_peff = P - flux_qof;
            flux_ea   = evap_3(st,S1,suzmax,Ep,delta_t);
            flux_qex  = saturation_1(flux_peff,S1,suzmax);
            flux_qv   = interflow_10(S1,kd,st*suzmax,suzmax-st*suzmax);
            flux_qb   = baseflow_4(q0,f,S2);

            % stores ODEs
            dS1 = flux_peff - flux_ea - flux_qex - flux_qv;
            dS2 = flux_qb - flux_qv;                                       % S2 is a deficit store
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_qof,  flux_peff, flux_ea,...
                      flux_qex, flux_qv, flux_qb];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end