classdef m_17_penman_4p_3s < MARRMoT_model
    % Class for pennman model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_17_penman_4p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 7;                                             % number of model fluxes
            obj.numParams = 4; 

            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1, 2000;    % smax, Maximum soil moisture storage [mm]
                             0, 1;       % phi, Fraction of direct runoff [-]
                             0, 1;       % gam, Evaporation reduction in lower zone [-]
                             0, 1];      % k1, Runoff coefficient [d-1]
            
            obj.StoreNames = ["S1" "S2" "S3"];                             % Names for the stores
            obj.FluxNames  = ["ea", "qex", "u1", "q12", "et", "u2", "q"];  % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 5];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [7];                                       % Index or indices of fluxes to add to Streamflow
            obj.StoreSigns  = [1 -1 1];                                    % Signs to give to stores (-1 is a deficit store), only needed for water balance

        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            smax  = theta(1);     % Maximum soil moisture storage [mm]
            phi   = theta(2);     % Fraction of direct runoff [-]
            gam   = theta(3);     % Evaporation reduction in lower zone [-]
            k1    = theta(4);     % Runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ea   = evap_1(S1,Ep,delta_t);
            flux_qex  = saturation_1(P,S1,smax);
            flux_u1   = split_1(phi,flux_qex);
            flux_q12  = split_1(1-phi,flux_qex);
            flux_et   = evap_16(gam,S2,S1,0.01,Ep,delta_t);
            flux_u2   = saturation_9(flux_q12,S2,0.01);
            flux_q    = baseflow_1(k1,S3);

            % stores ODEs
            dS1 = P       - flux_ea - flux_qex;
            dS2 = flux_et + flux_u2 - flux_q12;    
            dS3 = flux_u1 + flux_u2 - flux_q;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ea,  flux_qex, flux_u1,...
                      flux_q12, flux_et,  flux_u2, flux_q];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end