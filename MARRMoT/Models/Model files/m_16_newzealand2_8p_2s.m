classdef m_16_newzealand2_8p_2s < MARRMoT_model
% Class for hydrologic conceptual model: New Zealand model v2

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Atkinson, S. E., Sivapalan, M., Woods, R. A., & Viney, N. R. (2003). 
% Dominant physical controls on hourly flow predictions and the role of 
% spatial variability: Mahurangi catchment, New Zealand. Advances in Water 
% Resources, 26(3), 219–235. http://doi.org/10.1016/S0309-1708(02)00183-5

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
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
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"eint", "qtf", "veg", "ebs",...
                              "qse",  "qss", "qbf", "qt"};                 % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 3 4];                                   % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 8;                                         % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            s1max   = theta(1);     % Maximum interception storage [mm] 
            s2max   = theta(2);     % Maximum soil moisture storage [mm] 
            d       = theta(8);     % Routing delay [d]
            
            % maximum store values
            obj.store_max = [s1max, s2max];
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(d,delta_t);
            
            obj.uhs = {uh};
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
            uhs = obj.uhs;
            uh = uhs{1};
            
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
            flux_qt   = route(flux_qse + flux_qss + flux_qbf, uh);

            % stores ODEs
            dS1 = P - flux_eint - flux_qtf;
            dS2 = flux_qtf - flux_veg - flux_ebs -...
                  flux_qse - flux_qss - flux_qbf;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_eint, flux_qtf, flux_veg, flux_ebs,...
                      flux_qse,  flux_qss, flux_qbf, flux_qt];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % input fluxes to the unit hydrographs 
            fluxes = obj.fluxes(obj.t,:);
            flux_qse = fluxes(5);
            flux_qss = fluxes(6);
            flux_qbf = fluxes(7);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh = update_uh(uh, flux_qse + flux_qss + flux_qbf);
            obj.uhs = {uh};
        end
    end
end