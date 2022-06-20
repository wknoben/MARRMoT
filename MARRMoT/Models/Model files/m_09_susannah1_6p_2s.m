classdef m_09_susannah1_6p_2s < MARRMoT_model
% Class for hydrologic conceptual model: Susannah Brook v1

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Son, K., & Sivapalan, M. (2007). Improving model structure and reducing 
% parameter uncertainty in conceptual water balance models through the use 
% of auxiliary data. Water Resources Research, 43(1). 
% https://doi.org/10.1029/2006WR005032

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_09_susannah1_6p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 7;                                             % number of model fluxes
            obj.numParams = 6;
            
            obj.JacobPattern  = [1 0;...
                                 1 1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 2000;    % Sb, Maximum soil moisture storage [mm]
                             0.05, 0.95;    % Sfc, Wilting point as fraction of sb [-]
                             0.05, 0.95;    % M, Fraction forest [-]
                             1   , 50;      % a, Runoff coefficient [d] (should be > 0)
                             0.2 , 1;       % b, Runoff coefficient [-] (should be > 0)
                             0   , 1];      % r, Runoff coefficient [d-1] 
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"Ebs", "Evg", "Qse", "Qss",...
                              "qr",  "Qb",  "qt"};                         % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 2];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 7;                                         % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            sb  = theta(1);     % Maximum soil moisture storage [mm]
            sfc = theta(2);     % Wiliting point as fraction of sb [-]
            m   = theta(3);     % Fraction forest [-]
            a   = theta(4);     % Runoff coefficient [d]  
            b   = theta(5);     % Runoff coefficient [-]
            r   = theta(6);     % Runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
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
            flux_ebs  = evap_5(m,S1,sb,Ep,delta_t);
            flux_eveg = evap_6(m,sfc,S1,sb,Ep,delta_t);
            flux_qse  = saturation_1(P,S1,sb);
            flux_qss  = interflow_7(S1,sb,sfc,a,b,delta_t);
            flux_qr   = baseflow_1(r,flux_qss);
            flux_qb   = baseflow_2(S2,a,b,delta_t);
            flux_qt   = flux_qse + (flux_qss-flux_qr) + flux_qb;

            % stores ODEs
            dS1 = P - flux_ebs - flux_eveg - flux_qse - flux_qss;
            dS2 = flux_qr -  flux_qb;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_ebs, flux_eveg, flux_qse, flux_qss,...
                      flux_qr,  flux_qb,   flux_qt];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end