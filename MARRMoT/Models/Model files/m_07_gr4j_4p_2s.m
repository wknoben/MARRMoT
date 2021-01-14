classdef m_07_gr4j_4p_2s < MARRMoT_model
    % class for gr4j model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_07_gr4j_4p_2s(delta_t, theta)            
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 12;                                            % number of model fluxes
            obj.numParams = 4;
            
            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 2000;    % x1 [mm]
                            -20  , 20;      % x2 [mm/d]
                             1   , 300;     % x3 [mm]
                             .5  , 15];     % x4 [d]
            
            obj.StoreNames = ["S1" "S2" "S3" "S4" "S5"];                   % Names for the stores
            obj.FluxNames  = ["pn", "en", "ef", "ps", "es", "perc",...
                              "q9", "q1", "fr", "fq", "qr", "qt"];         % Names for the fluxes
            
            obj.Flux_Ea_idx = [3 5];                                       % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = 12;                                          % Index or indices of fluxes to add to Streamflow
            
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
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            x1 = theta(1);     % Maximum soil moisture storage [mm]
            x3 = theta(3);     % Maximum routing store storage [mm]
            x4 = theta(4);     % Flow delay [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = [x1, x3];
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh_q9 = uh_1_half(x4,delta_t);
            uh_q1 = uh_2_full(2*x4,delta_t);
            
            obj.uhs        = {uh_q9, uh_q1};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        % gr4j as implemented here is subtantially different that the
        % original MARRMoT: there, S1 is solved, then q1 and q9 are routed,
        % then S2 is solved, sequentially. Here, S1 and S2 are solved at
        % the same time, the results therefore are different.
        % I have implemented it in this way so that I can keep it
        % consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that S2 is balanced at all
        % steps, which is not the case in the original MARRMoT version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            x1      = theta(1);     % Maximum soil moisture storage [mm]
            x2      = theta(2);     % Water exchange coefficient [mm/d]
            x3      = theta(3);     % Maximum routing store storage [mm]
            x4      = theta(4);     % Flow delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh_q9 = uhs{1}; stf_q9 = stf{1};
            uh_q1 = uhs{2}; stf_q1 = stf{2};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pn   = max(P-Ep,0);
            flux_en   = max(Ep-P,0);
            flux_ef   = P - flux_pn;
            flux_ps   = saturation_4(S1,x1,flux_pn);
            flux_es   = evap_11(S1,x1,flux_en);
            flux_perc = percolation_3(S1,x1);
            flux_q9   = uh_q9(1).*.9.*(flux_pn - flux_ps + flux_perc) + stf_q9(1);
            flux_q1   = uh_q1(1).*.1.*(flux_pn - flux_ps + flux_perc) + stf_q1(1);
            flux_fr   = recharge_2(3.5,S2,x3,x2);
            flux_fq   = flux_fr;
            flux_qr   = baseflow_3(S2,x3);
            flux_qt   = flux_qr + max(flux_q1 + flux_fq,0);

            % stores ODEs
            dS1 = flux_ps - flux_es - flux_perc;
            dS2 = flux_q9 + flux_fr - flux_qr;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_pn, flux_en, flux_ef, flux_ps, flux_es, flux_perc,...
                      flux_q9, flux_q1, flux_fr, flux_fq, flux_qr, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh_q9 = uhs{1}; stf_q9 = stf{1};
            uh_q1 = uhs{2}; stf_q1 = stf{2};
            
            % input fluxes to the unit hydrographs  
            flux_pn   = fluxes(1);
            flux_ps   = fluxes(4);
            flux_perc = fluxes(6);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf_q9      = (uh_q9 .* (.9*(flux_pn - flux_ps + flux_perc))) + stf_q9;
            stf_q9      = circshift(stf_q9,-1);
            stf_q9(end) = 0;
            
            stf_q1      = (uh_q1 .* (.1*(flux_pn - flux_ps + flux_perc))) + stf_q1;
            stf_q1      = circshift(stf_q1,-1);
            stf_q1(end) = 0;
            
            obj.fluxes_stf = {stf_q9, stf_q1};
        end
    end
end