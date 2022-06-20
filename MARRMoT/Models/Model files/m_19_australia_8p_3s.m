classdef m_19_australia_8p_3s < MARRMoT_model
% Class for hydrologic conceptual model: Australia model

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Farmer, D., Sivapalan, M., & Jothityangkoon, C. (2003). Climate, soil, 
% and vegetation controls upon the variability of water balance in 
% temperate and semiarid landscapes: Downward approach to water balance 
% analysis. Water Resources Research, 39(2). 
% http://doi.org/10.1029/2001WR000328

    properties
        % model-specific attributes
    end
    methods
        
        % creator method

        function obj = m_19_australia_8p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 8;                                             % number of model fluxes
            obj.numParams = 8; 

            obj.JacobPattern  = [1,1,0;
                                 1,1,0;
                                 0,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 2000;    % Sb, Maximum soil moisture storage [mm]
                             0.05, 0.95;    % phi, Porosity [-]
                             0.01, 1.00;    % Sfc, Wilting point as fraction of sb [-]
                             0   , 1.00;    % alpha_ss, Subsurface flow constant [1/d]
                             1   , 5;       % beta_ss, Subsurface non-linearity constant [-]
                             0   , 1.00;    % k_deep, Groundwater recharge constant [d-1]
                             0   , 1.00;    % alpha_bf, Groundwater flow constant [d-1]
                             1   , 5];      % beta_bf, Groundwater non-linearity constant [-]
            
            obj.StoreNames = {"S1", "S2" "S3"};                             % Names for the stores
            obj.FluxNames  = {"eus", "rg",  "se", "esat",...
                              "qse", "qss", "qr", "qbf"};                  % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 4];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [5 6 8];                                   % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            sb       = theta(1);     % Maximum soil moisture storage [mm]
            phi      = theta(2);     % Porosity [-]
            fc       = theta(3);     % Field capacity as fraction of sb [-]
            alpha_ss = theta(4);     % Subsurface flow constant [1/d]
            beta_ss  = theta(5);     % Subsurface non-linearity constant [-]
            k_deep   = theta(6);     % Groundwater recharge constant [d-1]
            alpha_bf = theta(7);     % Groundwater flow constant [d-1]
            beta_bf  = theta(8);     % Groundwater non-linearity constant [-]
            
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
            flux_eus  = evap_7(S1,sb,Ep,delta_t);
            flux_rg   = saturation_1(P,S1,(sb-S2)*fc/phi);
            flux_se   = excess_1(S1,(sb-S2)*fc/phi,delta_t);
            flux_esat = evap_7(S2,sb,Ep,delta_t);
            flux_qse  = saturation_1(flux_rg+flux_se,S2,sb);
            flux_qss  = interflow_3(alpha_ss,beta_ss,S2,delta_t);
            flux_qr   = recharge_3(k_deep,S2);
            flux_qbf  = interflow_3(alpha_bf,beta_bf,S3,delta_t);

            % stores ODEs
            dS1 = P - flux_eus - flux_rg - flux_se;
            dS2 = flux_rg + flux_se - flux_esat - flux_qse - flux_qss - flux_qr;
            dS3 = flux_qr - flux_qbf;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_eus, flux_rg,  flux_se, flux_esat,...
                      flux_qse, flux_qss, flux_qr, flux_qbf];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end