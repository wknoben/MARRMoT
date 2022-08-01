classdef m_28_xinanjiang_12p_4s < MARRMoT_model
% Class for hydrologic conceptual model: Xinanjiang

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Jayawardena, A. W., & Zhou, M. C. (2000). A modified spatial soil moisture
% storage capacity distribution curve for the Xinanjiang model. Journal of 
% Hydrology, 227(1-4), 93–113. http://doi.org/10.1016/S0022-1694(99)00173-0
% 
% Zhao, R.-J. (1992). The Xinanjiang model applied in China. Journal of 
% Hydrology, 135(1-4), 371–381. http://doi.org/10.1016/0022-1694(92)90096-E

    properties
        % model-specific attributes
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_28_xinanjiang_12p_4s()
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 10;                                            % number of model fluxes
            obj.numParams = 12;
            
            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 0,1,1,0;
                                 0,1,0,1];                                 % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 1;           % aim,  Fraction impervious area [-]
                                -0.49, 0.49;    % a,    Tension water distribution inflection parameter [-]
                                0, 10;          % b,    Tension water distribution shape parameter [-]
                                1, 2000;        % stot, Total soil moisture storage (W+S) [mm]
                                0.01,0.99;      % fwm,  Fraction of Stot that is Wmax [-]
                                0.01,0.99;      % flm,  Fraction of wmax that is LM [-]
                                0.01, 0.99;     % c,    Fraction of LM for second evaporation change [-]
                                0, 10;          % ex,   Free water distribution shape parameter [-]
                                0, 1;           % ki,   Free water interflow parameter [d-1]
                                0, 1;           % kg,   Free water groundwater parameter [d-1]
                                0, 1;           % ci,   Interflow time coefficient [d-1]
                                0, 1];          % cg,   Baseflow time coefficient [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4"};                        % Names for the stores
            obj.FluxNames  = {"rb", "pi", "e",  "r",  "rs",...
                              "ri", "rg", "qs", "qi", "qg"};               % Names for the fluxes
            
            obj.FluxGroups.Ea = 3;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9 10];                                  % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta = obj.theta;
            stot = theta(4);     % Total soil moisture storage (W+S) [mm]
            fwmx = theta(5);     % Fraction of Stot that is Wmax [-]
            flm  = theta(6);     % Fraction of wmax that is LM [-]
            
            % auxiliary parameters
            wmax = fwmx*stot;    % Maximum tension water depth [mm]
            smax = (1-fwmx)*stot;% Maximum free water depth [mm]
            lm   = flm*wmax;     % Tension water threshold for evaporation change [mm]
            obj.aux_theta = [wmax, smax, lm];
            
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            aim  = theta(1);     % Fraction impervious area [-]
            a    = theta(2);     % Tension water distribution inflection parameter [-]
            b    = theta(3);     % Tension water distribution shape parameter [-]
            c    = theta(7);     % Fraction of LM for second evaporation change [-]
            ex   = theta(8);     % Free water distribution shape parameter [-]
            ki   = theta(9);     % Free water interflow parameter [d-1]
            kg   = theta(10);    % Free water baseflow parameter [d-1]
            ci   = theta(11);    % Interflow time coefficient [d-1]
            cg   = theta(12);    % Baseflow time coefficient [d-1]
            
            % auxiliary parameters
            aux_theta = obj.aux_theta;
            wmax = aux_theta(1); % Maximum tension water depth [mm]
            smax = aux_theta(2); % Maximum free water depth [mm]
            lm   = aux_theta(3); % Tension water threshold for evaporation change [mm]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_rb = split_1(aim,P);
            flux_pi = split_1(1-aim,P);
            flux_e  = evap_21(lm,c,S1,Ep,delta_t);
            flux_r  = saturation_14(a,b,S1,wmax,flux_pi);
            flux_rs = saturation_2(S2,smax,ex,flux_r);
            flux_ri = saturation_2(S2,smax,ex,S2*ki);
            flux_rg = saturation_2(S2,smax,ex,S2*kg);
            flux_qs = flux_rb + flux_rs;
            flux_qi = interflow_5(ci,S3);
            flux_qg = baseflow_1(cg,S4);
            
            % stores ODEs
            dS1 = flux_pi - flux_e  - flux_r;
            dS2 = flux_r  - flux_rs - flux_ri - flux_rg;
            dS3 = flux_ri - flux_qi;
            dS4 = flux_rg - flux_qg;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_rb flux_pi flux_e  flux_r  flux_rs ...
                      flux_ri flux_rg flux_qs flux_qi flux_qg];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end