classdef m_04_newzealand1_6p_1s < MARRMoT_model
    % Class for newzealand1 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_04_newzealand1_6p_1s()
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 5;                                             % number of model fluxes
            obj.numParams = 6;
            
            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [1   , 2000;   % Smax, Maximum soil moisture storage [mm]
                             0.05, 0.95;   % sfc, Field capacity as fraction of maximum soil moisture [-]
                             0.05, 0.95 ;  % m, Fraction forest [-]
                             0   , 1;      % a, Subsurface runoff coefficient [d-1]
                             1   , 5;      % b, Runoff non-linearity [-]
                             0   , 1];     % tcbf, Baseflow runoff coefficient [d-1]
            
            obj.StoreNames = ["S1"];                                       % Names for the stores
            obj.FluxNames  = ["veg", "ebs", "qse", "qss", "qbf"];          % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 2];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [3 4 5];                                   % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            s1max   = theta(1);     % Maximum soil moisture storage [mm] 
            sfc     = theta(2);     % Field capacity as fraction of maximum soil moisture [-]
            m       = theta(3);     % Fraction forest [-]
            a       = theta(4);     % Subsurface runoff coefficient [d-1]
            b       = theta(5);     % Runoff non-linearity [-]
            tcbf    = theta(6);     % Baseflow runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_veg  = evap_6(m,sfc,S1,s1max,Ep,delta_t);
            flux_ebs  = evap_5(m,S1,s1max,Ep,delta_t);
            flux_qse  = saturation_1(P,S1,s1max);
            flux_qss  = interflow_9(S1,a,sfc*s1max,b,delta_t);
            flux_qbf  = baseflow_1(tcbf,S1);

            % stores ODEs
            dS1 = P -flux_veg -flux_ebs -flux_qse -flux_qss -flux_qbf;
            
            % outputs
            dS = [dS1];
            fluxes = [flux_veg,  flux_ebs,...
                      flux_qse, flux_qss, flux_qbf];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end