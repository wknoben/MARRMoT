classdef m_29_hymod_5p_5s < MARRMoT_model
% Class for hydrologic conceptual model: HyMOD

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, Hoshin, 
% V., & Sorooshian, S. (2001). A framework for development and application 
% of hydrological models. Hydrology and Earth System Sciences, 5, 13–26.

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_29_hymod_5p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 8;                                             % number of model fluxes
            obj.numParams = 5;
            
            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 0,1,1,0,0;
                                 0,0,1,1,0;
                                 1,0,0,0,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1  , 2000;     % Smax, Maximum soil moisture storage [mm], 
                             0  , 10;       % b, Soil depth distribution parameter [-]
                             0  , 1;        % a, Runoff distribution fraction [-]
                             0  , 1;        % kf, fast flow time parameter [d-1]
                             0  , 1];       % ks, base flow time parameter [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"ea",  "pe",  "pf",  "ps",...
                              "qf1", "qf2", "qf3", "qs"};                  % Names for the fluxes
            
            obj.FluxGroups.Ea = 1;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [7 8];                                     % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)          
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            smax   = theta(1);     % Maximum soil moisture storage     [mm], 
            b      = theta(2);     % Soil depth distribution parameter [-]
            a      = theta(3);     % Runoff distribution fraction [-]
            kf     = theta(4);     % Fast runoff coefficient [d-1]
            ks     = theta(5);     % Slow runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ea  = evap_7(S1,smax,Ep,delta_t);
            flux_pe  = saturation_2(S1,smax,b,P);
            flux_pf  = split_1(a,flux_pe);
            flux_ps  = split_1(1-a,flux_pe);
            flux_qf1 = baseflow_1(kf,S2);
            flux_qf2 = baseflow_1(kf,S3);
            flux_qf3 = baseflow_1(kf,S4);
            flux_qs  = baseflow_1(ks,S5);

            % stores ODEs
            dS1 = P - flux_ea - flux_pe;
            dS2 = flux_pf - flux_qf1;
            dS3 = flux_qf1 - flux_qf2;
            dS4 = flux_qf2 - flux_qf3;
            dS5 = flux_ps - flux_qs;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_ea  flux_pe  flux_pf  flux_ps ...
                      flux_qf1 flux_qf2 flux_qf3 flux_qs];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end