classdef m_43_gsmsocont_12p_6s < MARRMoT_model
    % Class for gsmsocont model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_43_gsmsocont_12p_6s(delta_t, theta)
            obj.numStores = 6;                                             % number of model stores
            obj.numFluxes = 19;                                            % number of model fluxes
            obj.numParams = 12;

            obj.JacobPattern  = [1,0,0,0,0,0;
                                 1,1,0,0,0,0;
                                 1,0,1,0,0,0;
                                 0,0,0,1,0,0;
                                 0,0,0,1,1,0;
                                 0,0,0,1,1,1];                             % Jacobian matrix of model store ODEs
            
            obj.parRanges = [    0, 1;       % fice,   Fraction of catchment covered by glacier [-]
                                -3, 5;       % t0,     Threshold temperature for snowfall [oC]
                                 0, 20;      % asnow,  Degree-day factor for snow melt [mm/oC/d]
                                -3, 3;       % tm,     Threshold temperature for snow melt [oC]
                                 0, 1;       % ks,     Runoff coeficient for snow melt on glacier [d-1]
                                 0, 20;      % aice,   Degree-day factor for ice melt [mm/oC/d]
                                 0, 1;       % ki,     Runoff coeficient for ice melt on glacier [d-1]
                                 1, 2000;    % a,      Maximum soil moisture storage [mm]
                                 0, 10;      % x,      Evaporation non-linearity [-]
                                 0, 5;       % y,      Infiltration non-linearity [-]
                                 0, 1;       % ksl,    Runoff coefficient for baseflow [d-1]
                                 0, 1];      % beta,   Runoff coefficient for quick flow [mm^(4/3)/d]
        
            obj.StoreNames = ["S1" "S2" "S3" "S4" "S5" "S6"];              % Names for the stores
            obj.FluxNames  = ["pice", "pices", "picer", "mis", "pirs",...
                              "piri", "mii",   "qis",   "qii", "pni",...
                              "pnis", "pnir",  "mnis",  "peq", "peff",...
                              "pinf", "et",    "qsl",   "qqu"];            % Names for the fluxes
            
            obj.FluxGroups.Ea = [17];                                      % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9 18 19];                               % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.GlacierMelt = -7;                               % Index of flows from glacier melt (neg sign because it is input, differently from Q and Ea)
            
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
            fice  = theta(1);     % Fraction of catchment covered by glacier [-]
            t0    = theta(2);     % Threshold temperature for snowfall [oC]
            asnow = theta(3);     % Degree-day factor for snow melt [mm/oC/d]
            tm    = theta(4);     % Threshold temperature for snow melt [oC]
            ks    = theta(5);     % Runoff coeficient for snow melt on glacier [d-1]
            aice  = theta(6);     % Threshold temperature for ice melt [oC]
            ki    = theta(7);     % Runoff coeficient for ice melt on glacier [d-1]
            a     = theta(8);     % Maximum soil moisture storage [mm]
            x     = theta(9);     % Evaporation non-linearity [-]
            y     = theta(10);    % Infiltration non-linearity [-]
            ksl   = theta(11);    % Runoff coefficient for baseflow [d-1]
            beta  = theta(12);    % Runoff coefficient for quick flow [mm^(4/3)/d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            S6 = S(6);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pice    = split_1(fice,P);
            flux_pices   = snowfall_1(flux_pice,T,t0);
            flux_picer   = rainfall_1(flux_pice,T,t0);
            flux_mis     = melt_1(asnow,tm,T,S1,delta_t);
            flux_pirs    = saturation_9(flux_picer,S1,0.01);
            flux_piri    = effective_1(flux_picer,flux_pirs);
            flux_mii     = melt_3(aice,tm,T,Inf,S1,0.01,delta_t);
            flux_qis     = baseflow_1(ks,S2);
            flux_qii     = baseflow_1(ki,S3);
            flux_pni     = split_1(1-fice,P);
            flux_pnis    = snowfall_1(flux_pni,T,t0);
            flux_pnir    = rainfall_1(flux_pni,T,t0);
            flux_mnis    = melt_1(asnow,tm,T,S4,delta_t);
            flux_peq     = flux_pnir + flux_mnis;
            flux_peff    = infiltration_6(1,y,S5,a,flux_peq);
            flux_pinf    = effective_1(flux_peq,flux_peff);
            flux_et      = evap_19(1,x,S5,a,Ep,delta_t);
            flux_qsl     = baseflow_1(ksl,S5);
            flux_qqu     = interflow_3(beta,5/3,S6,delta_t);

            % stores ODEs
            dS1 = flux_pices - flux_mis;
            dS2 = flux_pirs  + flux_mis - flux_qis;
            dS3 = flux_piri  + flux_mii - flux_qii;
            dS4 = flux_pnis  - flux_mnis;
            dS5 = flux_pinf  - flux_et - flux_qsl;  
            dS6 = flux_peff  - flux_qqu; 
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6];
            fluxes = [flux_pice, flux_pices, flux_picer, flux_mis, flux_pirs,...
                      flux_piri, flux_mii,   flux_qis,   flux_qii, flux_pni,...
                      flux_pnis, flux_pnir,  flux_mnis,  flux_peq, flux_peff,...
                      flux_pinf, flux_et,    flux_qsl,   flux_qqu];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end