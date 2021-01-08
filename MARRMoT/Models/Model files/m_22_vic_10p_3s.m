classdef m_22_vic_10p_3s < MARRMoT_model
    % Class for vic model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
        aux_theta        % Auxiliary parameters
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_22_vic_10p_3s(delta_t, theta)
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 10;
            
            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0.1 , 5;        % ibar,     Mean interception capacity [mm]
                                0   , 1;        % idelta,   Seasonal interception change as fraction of mean [-]
                                1   , 365;      % ishift,   Maximum interception peak timing [-]
                                1   , 2000;     % stot,    Maximum soil moisture capacity [mm]
                                0.01, 0.99;     % fsm, Fraction of stot that constitutes maximum soil moisture smmax [-]
                                0   ,10;        % b,        Infiltration excess shape parameter [-]
                                0   , 1;        % k1,       Percolation time parameter [d-1]
                                0   ,10;        % c1,       Percolation non-linearity parameter [-]
                                0   ,1;         % k2,       Baseflow time parameter [d-1]
                                1   , 5];       % c2,       Baseflow non-linearity parameter
            
            obj.StoreNames = ["S1" "S2" "S3"];                             % Names for the stores
            obj.FluxNames  = ["ei",   "peff", "iex", "qie",  "inf", "et1",...
                              "qex1", "pc",   "et2", "qex2", "qb"];        % Names for the fluxes
            
            obj.Flux_Ea_idx = [1 6 9];                                     % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = [4 10 11];                                   % Index or indices of fluxes to add to Streamflow
            
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
            %Parameters
            theta = obj.theta;
            stot    = theta(4);     % Total available storage [mm]
            fsm     = theta(5);     % Fraction of stot that constitutes maximum soil mositure storage [-]
            
            % Auxiliary parameter
            smmax   = fsm*stot;     % Maximum soil moisture capacity [mm]
            gwmax   = (1-fsm)*stot; % Maximum groundwater storage [mm]
            tmax    = 365.25;       % Length of one growing cycle [d]
            obj.aux_theta = [smmax; gwmax; tmax];
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            ibar    = theta(1);     % Mean interception capacity [mm]
            idelta  = theta(2);     % Seasonal interception change as fraction of mean [-]
            ishift  = theta(3);     % Maximum interception peak timing [-]
            b       = theta(6);     % Infiltration excess shape parameter [-]
            k1      = theta(7);     % Percolation time parameter [d-1]
            c1      = theta(8);     % Percolation non-linearity parameter [-]
            k2      = theta(9);     % Baseflow time parameter [d-1]
            c2      = theta(10);    % Baseflow non-linearity parameter
            
            aux_theta = obj.aux_theta;
            smmax = aux_theta(1);
            gwmax = aux_theta(2);
            tmax = aux_theta(3);
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            aux_imax  = phenology_2(ibar,idelta,ishift,obj.t,tmax,delta_t);
            flux_ei   = evap_7(S1,aux_imax,Ep,delta_t);
            flux_peff = interception_1(P,S1,aux_imax);
            flux_iex  = excess_1(S1,aux_imax,delta_t);
            flux_qie  = saturation_2(S2,smmax,b,flux_peff+flux_iex);
            flux_inf  = effective_1(flux_peff+flux_iex,flux_qie);
            flux_et1  = evap_7(S2,smmax,max(0,Ep-flux_ei),delta_t);
            flux_qex1 = saturation_1(flux_inf,S2,smmax);
            flux_pc   = percolation_5(k1,c1,S2,smmax,delta_t);
            flux_et2  = evap_7(S3,gwmax,max(0,Ep-flux_ei-flux_et1),delta_t);
            flux_qex2 = saturation_1(flux_pc,S3,gwmax);
            flux_qb   = baseflow_5(k2,c2,S3,gwmax,delta_t);

            % stores ODEs
            dS1 = P - flux_ei - flux_peff - flux_iex;
            dS2 = flux_inf - flux_et1 - flux_qex1 - flux_pc;
            dS3 = flux_pc  - flux_et2 - flux_qex2 - flux_qb;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ei  flux_peff flux_iex  flux_qie ...
                      flux_inf flux_et1  flux_qex1 flux_pc ...
                      flux_et2 flux_qex2 flux_qb];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end