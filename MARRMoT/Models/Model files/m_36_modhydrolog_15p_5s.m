classdef m_36_modhydrolog_15p_5s < MARRMoT_model
% Class for hydrologic conceptual model: MODHYDROLOG

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Chiew, F. H. S. (1990). Estimating groundwater recharge using an
% integrated surface and groundwater model. University of Melbourne.
%
% Chiew, F., & McMahon, T. (1994). Application of the daily rainfall-runoff
% model MODHYDROLOG to 28 Australian catchments. Journal of Hydrology, 
% 153(1–4), 383–416. https://doi.org/10.1016/0022-1694(94)90200-3

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_36_modhydrolog_15p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 16;                                            % number of model fluxes
            obj.numParams = 15; 

            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 1,1,1,0,0;
                                 1,1,1,1,0;
                                 1,1,1,1,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0,     5;      % INSC, Maximum interception capacity, [mm]
                             0,     600;    % COEFF, Maximum infiltration loss parameter, [mm]
                             0,     15;     % SQ, Infiltration loss exponent, [-]
                             1,     2000;   % SMSC, Maximum soil moisture capacity, [mm]
                             0,     1;      % SUB, Proportionality constant, [-]
                             0,     1;      % CRAK, Proportionality constant, [-]
                             0,     20;     % EM, maximum plant-controled evap rate, [mm/d]
                             0,     50;     % DSC, Maximum depression capacity, [mm]
                             0,     1;      % ADS, Land fraction functioning as depression storage, [-]
                             0.99,  1;      % MD, Depression storage parameter, [-]
                             0,     0.5;    % VCOND, Leakage coefficient, [mm/d]
                           -10,     10;     % DLEV, Datum around which groundwater fluctuates relative to river bed, [mm]
                             0,     1;      % K1, Flow exchange parameter, [d-1] 
                             0,     1;      % K2, Flow exchange parameter, [d-1] 
                             0,     100];   % K3, Flow exchange parameter, [d-1]                  
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"Ei",   "EXC", "INF",  "INT",...
                              "REC",  "SMF", "Et",   "GWF",...
                              "TRAP", "Ed",  "DINF", "SEEP",...
                              "FLOW", "Q",   "RUN",  "SRUN"};              % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 7 10];                                  % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [14];                                      % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.GWseepage = 12;                                 % Index or GW seepage flux

        end
        
        % INITialisation function
        function obj = init(obj)
            obj.store_min(4) = -Inf;
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            insc    = theta(1);     % Maximum interception capacity, [mm]
            coeff   = theta(2);     % Maximum infiltration loss parameter, [-]
            sq      = theta(3);     % Infiltration loss exponent, [-]
            smsc    = theta(4);     % Maximum soil moisture capacity, [mm]
            sub     = theta(5);     % Proportionality constant, [-]
            crak    = theta(6);     % Proportionality constant, [-]
            em      = theta(7);     % Plant-controled maximum evaporation rate [mm/d]
            dsc     = theta(8);     % Maximum depression capacity, [mm]
            ads     = theta(9);     % Land fraction functioning as depression storage, [-]
            md      = theta(10);    % Depression storage parameter, [-], default = 1
            vcond   = theta(11);    % Leakage coefficient, [mm/d]
            dlev    = theta(12);    % Datum around which groundwater fluctuates relative to river bed, [mm]
            k1      = theta(13);    % Flow exchange parameter, [d-1] 
            k2      = theta(14);    % Flow exchange parameter, [d-1] 
            k3      = theta(15);    % Flow exchange parameter, [d-1] 
            
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
            flux_Ei     = evap_1(S1,Ep,delta_t);         
            flux_EXC    = interception_1(P,S1,insc);
            flux_INF    = infiltration_1(coeff,sq,S2,smsc,flux_EXC);
            flux_INT    = interflow_1(sub,S2,smsc,flux_INF);
            flux_REC    = recharge_1(crak,S2,smsc,flux_INF-flux_INT);
            flux_SMF    = flux_INF - flux_INT - flux_REC;
            flux_Et     = evap_2(em,S2,smsc,Ep,delta_t);
            flux_GWF    = saturation_1(flux_SMF,S2,smsc);
            flux_RUN    = flux_EXC - flux_INF;
            flux_TRAP   = depression_1(ads,md,S3,dsc,flux_RUN,delta_t);
            flux_Ed     = evap_1(S3,ads*Ep,delta_t);
            flux_DINF   = ads .* infiltration_2(coeff,sq,S2,smsc,flux_SMF,S3,delta_t);
            flux_SEEP   = exchange_3(vcond,S4,dlev);
            flux_SRUN   = flux_RUN - flux_TRAP;
            flux_FLOW   = exchange_1(k1,k2,k3,S4,flux_SRUN,delta_t);
            flux_Q      = baseflow_1(1,S5);

            % stores ODEs
            dS1 = P          - flux_Ei    - flux_EXC;
            dS2 = flux_SMF   + flux_DINF  - flux_Et    - flux_GWF;
            dS3 = flux_TRAP  - flux_Ed    - flux_DINF;
            dS4 = flux_REC   + flux_GWF   - flux_SEEP  - flux_FLOW;
            dS5 = flux_SRUN  + flux_INT   + flux_FLOW  - flux_Q; 
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_Ei,   flux_EXC, flux_INF,  flux_INT,...
                      flux_REC,  flux_SMF, flux_Et,   flux_GWF,...
                      flux_TRAP, flux_Ed,  flux_DINF, flux_SEEP,...
                      flux_FLOW, flux_Q,   flux_RUN,  flux_SRUN];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end