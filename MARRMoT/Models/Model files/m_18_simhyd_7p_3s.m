classdef m_18_simhyd_7p_3s < MARRMoT_model
% Class for hydrologic conceptual model: SimHyd

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Chiew, F. H. S., Peel, M. C., & Western, A. W. (2002). Application and 
% testing of the simple rainfall-runoff model SIMHYD. In V. P. Singh & D. 
% K. Frevert (Eds.), Mathematical Models of Small Watershed Hydrology (pp. 
% 335–367). Chelsea, Michigan, USA: Water Resources Publications LLC, USA.

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_18_simhyd_7p_3s()          
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 10;                                            % number of model fluxes
            obj.numParams = 7;                                             % number of model parameters
            
            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0,     5;      % INSC, Maximum interception capacity, [mm]
                             0,     600;    % COEFF, Maximum infiltration loss parameter, [mm]
                             0,     15;     % SQ, Infiltration loss exponent, [-]
                             1,     2000;   % SMSC, Maximum soil moisture capacity, [mm]
                             0,     1;      % SUB, Proportionality constant, [-]
                             0,     1;      % CRAK, Proportionality constant, [-]
                             0,     1];     % K, Slow flow time scale, [d-1]  
            
            obj.StoreNames = {"S1", "S2", "S3"};                           % Names for the stores
            obj.FluxNames  = {"Ei", "EXC", "INF", "INT", "REC",...
                              "Et", "GWF", "BAS", "SRUN", "Qt"};           % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 6];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 10;                                        % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function       
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            insc    = theta(1);     % Maximum interception capacity, [mm]
            coeff   = theta(2);     % Maximum infiltration loss parameter, [-]
            sq      = theta(3);     % Infiltration loss exponent, [-]
            smsc    = theta(4);     % Maximum soil moisture capacity, [mm]
            sub     = theta(5);     % Proportionality constant, [-],
            crak    = theta(6);     % Proportionality constant, [-]
            k       = theta(7);     % Slow flow time scale, [d-1]
            
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
            
            % fluxes functions
            % Original formulation using MARRMoT fluxes is slow on simhyd,
            % individual functions have been explicitly coded underneath.
            flux_Ei   = evap_1(S1,Ep,delta_t);
            flux_EXC  = interception_1(P,S1,insc);
            flux_INF  = infiltration_1(coeff,sq,S2,smsc,flux_EXC);
            flux_INT  = interflow_1(sub,S2,smsc,flux_INF);
            flux_REC  = recharge_1(crak,S2,smsc,(flux_INF-flux_INT));
            flux_Et   = evap_2(10,S2,smsc,Ep,delta_t);
            flux_GWF  = saturation_1((flux_INF-flux_INT-flux_REC),S2,smsc);
            flux_BAS  = baseflow_1(k,S3);
            flux_SRUN = flux_EXC - flux_INF; 
            flux_Qt   = flux_SRUN + flux_INT + flux_BAS;
            % flux_SMF is not reported in old MARRMoT, it is not reported
            % here either to keep results consistent
            flux_SMF  = flux_INF-flux_INT-flux_REC;
            
            % stores ODEs
            dS1  = P        - flux_Ei  - flux_EXC;
            dS2  = flux_SMF - flux_Et  - flux_GWF;   
            dS3  = flux_REC + flux_GWF - flux_BAS;
            
            % outputs
            dS = [dS1, dS2, dS3];
            fluxes = [flux_Ei, flux_EXC, flux_INF, flux_INT, flux_REC,...
                     flux_Et, flux_GWF, flux_BAS, flux_SRUN, flux_Qt];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end
