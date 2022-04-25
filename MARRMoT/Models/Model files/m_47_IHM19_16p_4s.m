classdef m_47_IHM19_16p_4s < MARRMoT_model
% Class for hydrologic conceptual model: Forellenbach model
    
% Copyright (C) 2021 Clara Brandes, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Brandes, C. (2020). Erstellung eines Konzeptionellen
% Hochwasserabfluss-modells für das Einzugsgebiet des Forellenbachs, NP
% Bayerischer Wald. MSc Thesis. Technische Univerisität Dresden, Germany.

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_47_IHM19_16p_4s()            
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 18;                                            % number of model fluxes
            obj.numParams = 16;
            
            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 1,1,1,0;
                                 0,1,1,1];                                 % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   2,     5;      % 01 SIMAX, maximum interception storage [mm]
                                0.9,   1;      % 02 A, splitting coeffcient for excess precipitation [-]
                                0.4,   0.95;   % 03 FF, forest fraction [-]
                                0.05,  5;      % 04 SMPMAX, maximum storage macropores [mm]
                                0,     1;      % 05 CQMP, runoff time parameter (fast/slow runnoff) first soil layer [1/d]
                                1,     5;      % 06 XQMP, runoff scale parameter first soil layer [-]     
                              400,   600;      % 07 SS1MAX, maximum soil moisture storage first soil layer [mm]
                                0.3,   0.7;    % 08 FCCS1, field capacity as fraction of maximum storage first soil layer [-]
                                0,  1000;      % 09 CFS1, maximum infiltration rate first soil layer [mm/d]
                                0,    15;      % 10 XFS1, infiltration loss exponent first soil layer [-]
                                0,     1;      % 11 CQS1, runoff time parameter for (fast/slow runnoff) first soil layer [1/d]
                                1,     5;      % 12 XQS1, runoff scale parameter first soil layer [-]
                              300,   500;      % 13 SS2MAX, maximum soil moisture storage second soil layer [mm]
                                0,     1;      % 14 CQS2, runoff time parameter for (fast/slow runnoff) second soil layer [1/d]
                                1,     5;      % 15 XQS2, runoff scale parameter second soil layer [-]
                                0.01,  5];     % 16 D, Flow delay before surface runoff [d]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4"};                        % Names for the stores
            obj.FluxNames  = {"ei",    "pex", "pexmp",  "pexs1", "fmp",...
                              "qexmp", "qmp", "pqexsl", "fs1",   "etas1",...
                              "qs1",   "q0",  "q0r",    "qmps1", "pc",...
                              "qh",    "qs2", "qgw"};                      % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 10];                                    % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [13 16 17 18];                             % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.GW = -18;                                       % Index of GW runoff flow.
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            SIMAX   = theta(1);                         % maximum interception storage [mm]
            SMPMAX  = theta(4);                         % maximum storage macropores [mm]
            SS1MAX  = theta(7);                         % maximum soil moisture storage first soil layer [mm]
            SS2MAX  = theta(13);                        % maximum soil moisture storage second soil layer [mm]
            D       = theta(16);                        % Flow delay before surface runoff [d]
            
            % initial store values
            obj.S0 = 0.8 * [SIMAX;                   % Initial interception storage
                            SMPMAX;                  % Initial macropore storage 
                            SS1MAX;                  % Initial soil moisture storage 1 
                            SS2MAX];                 % Initial soil moisture storage 2
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh_q0r = uh_4_full(D,delta_t);
            
            obj.uhs = {uh_q0r};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            SIMAX   = theta(1);                         % maximum interception storage [mm]
            A       = theta(2);                         % splitting coeffcient for excess precipitation [-]
            FF      = theta(3);                         % forest fraction [-]
            SMPMAX  = theta(4);                         % maximum storage macropores [mm]
            CQMP    = theta(5);                         % runoff time parameter (fast/slow runnoff) first soil layer [1/d]
            XQMP    = theta(6);                         % runoff scale parameter first soil layer [-]
            SS1MAX  = theta(7);                         % maximum soil moisture storage first soil layer [mm]
            FCCS1   = theta(8);                         % field capacity coefficient fist soil layer [-]
            CFS1    = theta(9);                         % maximum infiltration rate first soil layer [-]
            XFS1    = theta(10);                        % infiltration loss exponent first soil layer [-]
            CQS1    = theta(11);                        % runoff time parameter for (fast/slow runnoff) first soil layer [1/d]
            XQS1    = theta(12);                        % runoff scale parameter first soil layer [-]
            SS2MAX  = theta(13);                        % maximum soil moisture storage second soil layer [mm]
            CQS2    = theta(14);                        % runoff time parameter for (fast/slow runnoff) second soil layer [1/d]
            XQS2    = theta(15);                        % runoff scale parameter second soil layer [-]
            
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q0r = uhs{1};
            
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
            flux_ei      = evap_1(S1,Ep,delta_t);
            flux_pex     = interception_1(P,S1,SIMAX);
            flux_pexmp   = split_1(A,flux_pex);
            flux_pexs1   = split_2(A,flux_pex);
            flux_fmp     = infiltration_3(flux_pexmp,S2,SMPMAX);
            flux_qexmp   = flux_pexmp - flux_fmp;
            flux_qmp     = interflow_3(CQMP,XQMP,S2,delta_t);
            flux_pqexs1  = flux_pexs1 + flux_qexmp;
            flux_fs1     = infiltration_7(CFS1,XFS1,S3,SS1MAX,flux_pqexs1);
            flux_etas1   = evap_23(FF,FCCS1,S3,SS1MAX,Ep,delta_t);  
            flux_qs1     = interflow_9(S3,CQS1,FCCS1*SS1MAX,XQS1,delta_t);
            flux_q0      = flux_pqexs1 - flux_fs1;
            flux_q0r     = route(flux_q0, uh_q0r);
            flux_qmps1   = flux_qmp + flux_qs1;
            flux_pc      = infiltration_3(flux_qmps1,S4,SS2MAX);
            flux_qh      = flux_qmps1 - flux_pc;
            flux_qs2     = interflow_3(CQS2,XQS2,S4,delta_t);
            flux_qgw     = 0.0195;

            % stores ODEs
            dS1 = P - flux_ei - flux_pex;
            dS2 = flux_fmp - flux_qmp;
            dS3 = flux_fs1 - flux_etas1 - flux_qs1;
            dS4 = flux_pc - flux_qs2;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_ei, flux_pex, flux_pexmp, flux_pexs1,...
                      flux_fmp, flux_qexmp, flux_qmp, flux_pqexs1,...
                      flux_fs1, flux_etas1, flux_qs1, flux_q0,...
                      flux_q0r, flux_qmps1, flux_pc , flux_qh,...
                      flux_qs2, flux_qgw];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q0r = uhs{1};
            
            % input fluxes to the unit hydrographs
            fluxes = obj.fluxes(obj.t,:);
            flux_q0   = fluxes(12);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh_q0r = update_uh(uh_q0r, flux_q0);
            obj.uhs = {uh_q0r};
        end
    end
end