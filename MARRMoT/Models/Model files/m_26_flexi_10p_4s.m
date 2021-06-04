classdef m_26_flexi_10p_4s < MARRMoT_model
    % Class for flexi model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_26_flexi_10p_4s(delta_t, theta)
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 10; 

            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 1,1,1,0;
                                 1,1,0,1];                                 % Jacobian matrix of model store ODEs
                              
            obj.parRanges = [1,2000;        % URmax, Maximum soil moisture storage [mm]
                             0, 10;         % beta, Unsaturated zone shape parameter [-]
                             0, 1;          % D, Fast/slow runoff distribution parameter [-]
                             0, 20;         % PERCmax, Maximum percolation rate [mm/d]
                             0.05, 0.95;    % Lp, Wilting point as fraction of s1max [-]
                             1, 5;          % Nlagf, Flow delay before fast runoff [d]
                             1, 15;         % Nlags, Flow delay before slow runoff [d]
                             0, 1;          % Kf, Fast runoff coefficient [d-1]
                             0, 1;          % Ks, Slow runoff coefficient [d-1]
                             0, 5];         % Imax, Maximum interception storage [mm]
            
            obj.StoreNames = ["S1" "S2" "S3" "S4"];                        % Names for the stores
            obj.FluxNames  = ["peff", "ei", "ru",  "eur", "ps",...
                              "rf",   "rs", "rfl", "rsl", "qf", "qs"];     % Names for the fluxes
            
            obj.FluxGroups.Ea = [2 4];                                       % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [10 11];                                     % Index or indices of fluxes to add to Streamflow
            
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
            nlagf   = theta(6);     % Flow delay before fast runoff [d]
            nlags   = theta(7);     % Flow delay before slow runoff [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh_f = uh_3_half(nlagf,delta_t);
            uh_s = uh_3_half(nlags,delta_t);
            
            obj.uhs        = {uh_f, uh_s};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        % flexb as implemented here is subtantially different that the
        % original MARRMoT: there, S1 and S2 are solved, then rsl and rfl
        % are routed, then S3 and S4 are solved, sequentially. Here, all
        % stores are solved at the same time, the results therefore are
        % different. I have implemented it in this way so that I can keep
        % it consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that S3 and S4 are balanced at
        % all steps, which is not the case in the original MARRMoT version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            smax    = theta(1);     % Maximum soil moisture storage [mm]
            beta    = theta(2);     % Unsaturated zone shape parameter [-]
            d       = theta(3);     % Fast/slow runoff distribution parameter [-]
            percmax = theta(4);     % Maximum percolation rate [mm/d]
            lp      = theta(5);     % Wilting point as fraction of s1max [-]
            kf      = theta(8);     % Fast runoff coefficient [d-1]
            ks      = theta(9);     % Slow runoff coefficient [d-1]
            imax    = theta(10);    % Maximum interception storage [mm]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh_f = uhs{1}; stf_f = stf{1};
            uh_s = uhs{2}; stf_s = stf{2};
            
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
            flux_peff= interception_1(P,S1,imax);
            flux_ei  = evap_1(S1,Ep,delta_t);
            flux_ru  = saturation_3(S2,smax,beta,flux_peff);
            flux_eur = evap_3(lp,S2,smax,Ep,delta_t);
            flux_ps  = percolation_2(percmax,S2,smax,delta_t);
            flux_rf  = split_1(1-d,flux_peff-flux_ru);
            flux_rs  = split_1(d,flux_peff-flux_ru);
            flux_rfl = uh_f(1).*(flux_rf) + stf_f(1);
            flux_rsl = uh_s(1).*(flux_ps + flux_rs) + stf_s(1);
            flux_qf  = baseflow_1(kf,S3);
            flux_qs  = baseflow_1(ks,S4);

            % stores ODEs
            dS1 = P - flux_peff - flux_ei;
            dS2 = flux_ru - flux_eur - flux_ps;
            dS3 = flux_rfl - flux_qf;    
            dS4 = flux_rsl - flux_qs;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_peff, flux_ei, flux_ru,  flux_eur, flux_ps,...
                      flux_rf,   flux_rs, flux_rfl, flux_rsl, flux_qf, flux_qs];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh_f = uhs{1}; stf_f = stf{1};
            uh_s = uhs{2}; stf_s = stf{2};
            
            % input fluxes to the unit hydrographs  
            flux_ps = fluxes(5);
            flux_rf = fluxes(6);
            flux_rs = fluxes(7);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf_f      = (uh_f .* (flux_rf)) + stf_f;
            stf_f      = circshift(stf_f,-1);
            stf_f(end) = 0;
            
            stf_s      = (uh_s .* (flux_ps + flux_rs)) + stf_s;
            stf_s      = circshift(stf_s,-1);
            stf_s(end) = 0;
            
            obj.fluxes_stf = {stf_f, stf_s};
        end
    end
end