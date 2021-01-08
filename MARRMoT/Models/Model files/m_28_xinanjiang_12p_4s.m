classdef m_28_xinanjiang_12p_4s < MARRMoT_model
    % Class for xinanjiang model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_28_xinanjiang_12p_4s(delta_t, theta)
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 10;                                            % number of model fluxes
            obj.numParams = 12;
            
            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 0,1,1,0;
                                 0,0,1,1];                                 % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 1;           % aim,  Fraction impervious area [-]
                                -0.49, 0.49;    % a,    Tension water distribution inflection parameter [-]
                                0, 10;          % b,    Tension water distribution shape parameter [-]
                                1, 2000;        % stot, Total soil moisture storage (W+S) [mm]
                                0.01,0.99;      % fwm,  Fraction of Stot that is Wmax [-]
                                0.01,0.99;      % flm,  Fraction of wmax that is LM [-]
                                0.01, 0.99;     % c,    Fraction of LM for second evaporation change [-]
                                0, 10;          % ex,   Free water distribution shape parameter [-]
                                0, 1;           % ki,   Free water interflow parameter [d-1]
                                0, 1;           % kg,   Free water groundwater parameter [d-1]
                                0, 1;           % ci,   Interflow time coefficient [d-1]
                                0, 1];          % cg,   Baseflow time coefficient [d-1]
            
            obj.StoreNames = ["S1" "S2" "S3" "S4"];                        % Names for the stores
            obj.FluxNames  = ["rb", "pi", "e",  "r",  "rs",...
                              "ri", "rg", "qs", "qi", "qg"];               % Names for the fluxes
            
            obj.Flux_Ea_idx = 3;                                           % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = [8 9 10];                                    % Index or indices of fluxes to add to Streamflow
            
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
            theta = obj.theta;
            stot = theta(4);     % Total soil moisture storage (W+S) [mm]
            fwmx = theta(5);     % Fraction of Stot that is Wmax [-]
            flm  = theta(6);     % Fraction of wmax that is LM [-]
            
            % auxiliary parameters
            wmax = fwmx*stot;    % Maximum tension water depth [mm]
            smax = (1-fwmx)*stot;% Maximum free water depth [mm]
            lm   = flm*wmax;     % Tension water threshold for evaporation change [mm]
            obj.aux_theta = [wmax, smax, lm];
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            aim  = theta(1);     % Fraction impervious area [-]
            a    = theta(2);     % Tension water distribution inflection parameter [-]
            b    = theta(3);     % Tension water distribution shape parameter [-]
            c    = theta(7);     % Fraction of LM for second evaporation change [-]
            ex   = theta(8);     % Free water distribution shape parameter [-]
            ki   = theta(9);     % Free water interflow parameter [d-1]
            kg   = theta(10);    % Free water baseflow parameter [d-1]
            ci   = theta(11);    % Interflow time coefficient [d-1]
            cg   = theta(12);    % Baseflow time coefficient [d-1]
            
            % auxiliary parameters
            aux_theta = obj.aux_theta;
            wmax = aux_theta(1); % Maximum tension water depth [mm]
            smax = aux_theta(2); % Maximum free water depth [mm]
            lm   = aux_theta(3); % Tension water threshold for evaporation change [mm]
            
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
            flux_rb = split_1(aim,P);
            flux_pi = split_1(1-aim,P);
            flux_e  = evap_21(lm,c,S1,Ep,delta_t);
            flux_r  = saturation_14(a,b,S1,wmax,flux_pi);
            flux_rs = saturation_2(S2,smax,ex,flux_r);
            flux_ri = saturation_2(S2,smax,ex,S2*ki);
            flux_rg = saturation_2(S2,smax,ex,S2*kg);
            flux_qs = flux_rb + flux_rs;
            flux_qi = interflow_5(ci,S3);
            flux_qg = baseflow_1(cg,S4);
            
            % stores ODEs
            dS1 = flux_pi - flux_e  - flux_r;
            dS2 = flux_r  - flux_rs - flux_ri - flux_rg;
            dS3 = flux_ri - flux_qi;
            dS4 = flux_rg - flux_qg;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_rb flux_pi flux_e  flux_r  flux_rs ...
                      flux_ri flux_rg flux_qs flux_qi flux_qg];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end