classdef m_08_us1_5p_2s < MARRMoT_model
    % Class for us1 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_08_us1_5p_2s(delta_t, theta)
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 9;                                             % number of model fluxes
            obj.numParams = 5; 

            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0   , 1;       % Alpha_ei, Fraction of intercepted rainfall [-]
                             0.05, 0.95;    % M, Fraction forest [-]
                             1   , 2000;    % Smax, Maximum soil moisture [mm]
                             0.05, 0.95;    % fc, Field capacity as fraction of Smax [-]
                             0   , 1];      % Alpha_ss, Subsurface routing delay [d-1]
            
            obj.StoreNames = ["S1" "S2"];                                  % Names for the stores
            obj.FluxNames  = ["eusei",  "eusveg", "eusbs", "esatveg",...
                              "esatbs", "rg",     "se",    "qse", "qss"];  % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 2 3 4 5];                               % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9];                                     % Index or indices of fluxes to add to Streamflow
            
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
            theta = obj.theta;
            alpha_ei    = theta(1);     % Fraction of intercepted rainfall [-]
            m           = theta(2);     % Fraction forest [-]
            smax        = theta(3);     % Maximum soil moisture [mm]
            fc          = theta(4);     % Field capacity as fraction of smax [-]
            alpha_ss    = theta(5);     % Subsurface routing delay [d-1]
            
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
            flux_eusei   = interception_3(alpha_ei,P);
            flux_eusveg  = evap_8(S1,S2,m,fc*(smax-S2),Ep,delta_t);
            flux_eusbs   = evap_9(S1,S2,m,smax,Ep,delta_t);
            flux_esatveg = evap_10(m,S2,S1+S2,Ep,delta_t);
            flux_esatbs  = evap_5(m,S2,S1+S2,Ep,delta_t);
            flux_rg      = saturation_1(P,S1,fc*(smax-S2));
            flux_se      = excess_1(S1,fc*(smax-S2),delta_t);
            flux_qse     = saturation_1(flux_rg+flux_se,S2,smax);
            flux_qss     = baseflow_1(alpha_ss,S2);

            % stores ODEs
            dS1 =   P     - flux_eusei - flux_eusveg  - flux_eusbs  - flux_rg  - flux_se;
            dS2 = flux_rg + flux_se    - flux_esatveg - flux_esatbs - flux_qse - flux_qss;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_eusei,   flux_eusveg, flux_eusbs, ...
                      flux_esatveg, flux_esatbs, flux_rg,...
                      flux_se,      flux_qse,    flux_qss];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
        end
    end
end