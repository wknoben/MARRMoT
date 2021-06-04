classdef m_25_tcm_6p_4s < MARRMoT_model
    % Class for tcm model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_25_tcm_6p_4s(delta_t, theta)
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 6;
            
            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 1,1,1,0
                                 0,0,1,1];                                 % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 1;     % phi, Fraction preferential recharge [-]
                                1, 2000;  % rc, Maximum soil moisture depth [mm]
                                0, 1;     % gam, Fraction of Ep reduction with depth [-]
                                0, 1;     % k1, Runoff coefficient [d-1]
                                0, 1;     % fa, Fraction of mean(P) that forms abstraction rate [mm/d]
                                0, 1];    % k2, Runoff coefficient [mm-1 d-1]
            
            obj.StoreNames = ["S1" "S2" "S3" "S4"];                        % Names for the stores
            obj.FluxNames  = ["en",  "ea",   "et",   "pn",  "pby",...
                              "pin", "qex1", "qex2", "qux", "a", "q"];     % Names for the fluxes
            
            obj.FluxGroups.Ea_idx = [1 2 3];                               % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q_idx  = 11;                                    % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.Abstraction = 10;                               % Index or abstraction flux (just needed for water balance)
            obj.StoreSigns  = [1 -1 1 1];                                  % Signs to give to stores (-1 is a deficit store), only needed for water balance

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
            % Original MARRMoT has fa = theta(5) and ca = fa * mean(P) as
            % an auxiliary parameter. Here ca = theta(5) is a parameter by
            % itself as described in Moore & Bell (2001) in order to keep
            % it consisted with the other models and remove the dependence
            % on the entire record of input rainfall.
            theta   = obj.theta;
            phi   = theta(1);     % Fraction preferential recharge [-]
            rc    = theta(2);     % Maximum soil moisture depth [mm]
            gam   = theta(3);     % Fraction of Ep reduction with depth [-]
            k1    = theta(4);     % Runoff coefficient [d-1]
            ca    = theta(5);     % Abstraction rate [mm/d]
            k2    = theta(6);     % Runoff coefficient [mm-1 d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pn     = effective_1(P,Ep);
            flux_en     = P - flux_pn;
            flux_pby    = split_1(phi,flux_pn);
            flux_pin    = split_1(1-phi,flux_pn);
            flux_ea     = evap_1(S1,Ep,delta_t);
            flux_et     = evap_16(gam,S2,S1,0.01,Ep,delta_t);
            flux_qex1   = saturation_1(flux_pin,S1,rc);
            flux_qex2   = saturation_9(flux_qex1,S2,0.01);
            flux_quz    = baseflow_1(k1,S3);
            flux_a      = abstraction_1(ca);
            flux_q      = baseflow_6(k2,0,S4,delta_t);
            
            % stores ODEs
            dS1 = flux_pin  - flux_ea   - flux_qex1;
            dS2 = flux_et   + flux_qex2 - flux_qex1;    
            dS3 = flux_qex2 + flux_pby  - flux_quz;
            dS4 = flux_quz  - flux_a    - flux_q;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_pn, flux_en, flux_pby,  flux_pin,...
                      flux_ea, flux_et, flux_qex1, flux_qex2,...
                      flux_quz, flux_a, flux_q]; 
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end