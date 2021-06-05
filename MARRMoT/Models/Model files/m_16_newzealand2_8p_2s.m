classdef m_16_newzealand2_8p_2s < MARRMoT_model
    % Class for newzealand2 model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_16_newzealand2_8p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 8;                                             % number of model fluxes
            obj.numParams = 8;
            
            obj.JacobPattern  = [1,0;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0   , 5;      % Maximum interception storage [mm] 
                             1   , 2000;   % Smax, Maximum soil moisture storage [mm]
                             0.05, 0.95;   % sfc, Field capacity as fraction of maximum soil moisture [-]
                             0.05, 0.95 ;  % m, Fraction forest [-]
                             0   , 1;      % a, Subsurface runoff coefficient [d-1]
                             1   , 5;      % b, Runoff non-linearity [-]
                             0   , 1;      % tcbf, Baseflow runoff coefficient [d-1]
                             1   , 120];   % Routing time delay [d]
            
            obj.StoreNames = ["S1" "S2"];                                  % Names for the stores
            obj.FluxNames  = ["eint", "qtf", "veg", "ebs",...
                              "qse",  "qss", "qbf", "qt"];                 % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 3 4];                                   % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 8;                                         % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            d = theta(8);     % Routing delay [d]
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(d,delta_t);
            
            obj.uhs        = {uh};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            s1max   = theta(1);     % Maximum interception storage [mm] 
            s2max   = theta(2);     % Maximum soil moisture storage [mm] 
            sfc     = theta(3);     % Field capacity as fraction of maximum soil moisture [-]
            m       = theta(4);     % Fraction forest [-]
            a       = theta(5);     % Subsurface runoff coefficient [d-1]
            b       = theta(6);     % Runoff non-linearity [-]
            tcbf    = theta(7);     % Baseflow runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_eint = evap_1(S1,Ep,delta_t);
            flux_qtf  = interception_1(P,S1,s1max);
            flux_veg  = evap_6(m,sfc,S2,s2max,Ep,delta_t);
            flux_ebs  = evap_5(m,S2,s2max,Ep,delta_t);
            flux_qse  = saturation_1(flux_qtf,S2,s2max);
            flux_qss  = interflow_9(S2,a,sfc*s2max,b,delta_t);
            flux_qbf  = baseflow_1(tcbf,S2);
            flux_qt   = uh(1) * (flux_qse + flux_qss + flux_qbf) + stf(1);

            % stores ODEs
            dS1 = P - flux_eint - flux_qtf;
            dS2 = flux_qtf - flux_veg - flux_ebs -...
                  flux_qse - flux_qss - flux_qbf;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_eint, flux_qtf, flux_veg, flux_ebs,...
                      flux_qse,  flux_qss, flux_qbf, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % input fluxes to the unit hydrographs
            flux_qse = fluxes(5);
            flux_qss = fluxes(6);
            flux_qbf = fluxes(7);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf      = (uh .* (flux_qse + flux_qss + flux_qbf)) + stf;
            stf      = circshift(stf,-1);
            stf(end) = 0;
            
            obj.fluxes_stf = {stf};
        end
    end
end