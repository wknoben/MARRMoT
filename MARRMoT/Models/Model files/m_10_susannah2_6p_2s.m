classdef m_10_susannah2_6p_2s < MARRMoT_model
% Class for hydrologic conceptual model: Susannah Brook model v2

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
        function obj = m_10_susannah2_6p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 8;                                             % number of model fluxes
            obj.numParams = 6;
            
            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 2000;    % Sb, Maximum soil moisture storage [mm]
                             0.05, 0.95;    % phi, Porosity [-]
                             0.05, 0.95;    % Sfc, Wilting point as fraction of sb [-]
                             0   , 1 ;      % r, Fraction of runoff coefficient [-]
                             0   , 1 ;      % c, Runoff coefficient [d-1] (should be > 0)
                             1   , 5];      % d, Runoff coefficient [-] (should be > 0)
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"eus", "rg",  "se", "esat",...
                              "qse", "qss", "qr", "qt"};                   % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 4];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 8;                                         % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.GWsink = 7;                                     % Index or groundwater sink flux
            
        end
        
        % INITialisation function
        function obj = init(obj)            
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            sb    = theta(1);     % Maximum soil moisture storage [mm]
            phi   = theta(2);     % Porosity [-]
            fc    = theta(3);     % Field capacity as fraction of sb [-]
            r     = theta(4);     % Fraction of recharge coefficient [-]
            c     = theta(5);     % Subsurface flow constant [1/d]
            d     = theta(6);     % Subsurface flow constant [-]
            
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
            flux_eus  = evap_7(S1,sb,Ep,delta_t);
            flux_rg   = saturation_1(P,S1,(sb-S2)*fc/phi);
            flux_se   = excess_1(S1,(sb-S2)*fc/phi,delta_t);
            flux_esat = evap_7(S2,sb,Ep,delta_t);
            flux_qse  = saturation_1(flux_rg+flux_se,S2,sb);
            flux_qss  = interflow_3((1-r)*c,d,S2,delta_t);
            flux_qr   = interflow_3(r*c,d,S2,delta_t);
            flux_qt   = flux_qse + flux_qss;

            % stores ODEs
            dS1 = P - flux_eus - flux_rg - flux_se;
            dS2 = flux_rg + flux_se - flux_esat - ...
                  flux_qse - flux_qss - flux_qr;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_eus, flux_rg,  flux_se, flux_esat,...
                      flux_qse, flux_qss, flux_qr, flux_qt];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end