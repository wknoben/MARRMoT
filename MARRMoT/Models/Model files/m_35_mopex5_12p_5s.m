classdef m_35_mopex5_12p_5s < MARRMoT_model
% Class for hydrologic conceptual model: MOPEX-5

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Ye, S., Yaeger, M., Coopersmith, E., Cheng, L., & Sivapalan, M. (2012). 
% Exploring the physical controls of regional patterns of flow duration 
% curves - Part 2: Role of seasonality, the regime curve, and associated 
% process controls. Hydrology and Earth System Sciences, 16(11), 4447–4465.
% http://doi.org/10.5194/hess-16-4447-2012

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_35_mopex5_12p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 13;                                            % number of model fluxes
            obj.numParams = 12;
            
            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 0,1,1,0,0;
                                 0,1,1,1,0;
                                 0,0,1,0,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [-3  , 3;       % tcrit, Snowfall & snowmelt temperature [oC]
                             0   , 20;      % ddf, Degree-day factor for snowmelt [mm/oC/d]
                             1   , 2000;    % Sb1, Maximum soil moisture storage [mm]
                             0   , 1 ;      % tw, Groundwater leakage time [d-1]
                             0   , 1 ;      % I_alpha, Intercepted fraction of Pr [-]
                             1   , 365;     % I_s, Maximum Leaf Area Index timing [d]
                             -10 , 0;       % tmin, Growing Season Index minimum temperature [oC]
                             1   , 20;      % trange, Growing Season Index temperature range [oC]
                             0   , 1 ;      % tu, Slow flow routing response time [d-1]
                             0.05, 0.95;    % se, Root zone storage capacity as fraction of Sb2 [-]
                             1   , 2000;    % Sb2, Root zone storage capacity [mm]
                             0   , 1];      % tc, Mean residence time [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"epc", "ps", "pr",  "qn",  "et1", "i",...
                              "q1f", "qw", "et2", "q2f", "q2u", "qf", "qs"};% Names for the fluxes
            
            obj.FluxGroups.Ea = [5 6 9];                                   % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [12 13];                                   % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            tcrit   = theta(1);     % Snowfall & snowmelt temperature [oC]
            ddf     = theta(2);     % Degree-day factor for snowmelt [mm/oC/d]
            s2max   = theta(3);     % Maximum soil moisture storage [mm]
            tw      = theta(4);     % Groundwater leakage time [d-1]
            i_alpha = theta(5);     % Intercepted fraction of Pr [-]
            i_s     = theta(6);     % Maximum Leaf Area Index timing [d]
            tmin    = theta(7);     % Growing Season Index minimum temperature
            trange  = theta(8);     % Growing Season Index temperature range
            tu      = theta(9);     % Slow flow routing response time [d-1]
            se      = theta(10);    % Root zone storage capacity as fraction of s3max [-]
            s3max   = theta(11);    % Maximum groundwater storage [mm]
            tc      = theta(12);    % Mean residence time [d-1]
            
            % auxiliary parameters
            tmax   = 365.25;        % Duration of seasonal cycle [d]
            
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
            flux_epc  = phenology_1(T,tmin,tmin+trange,Ep);
            flux_ps   = snowfall_1(P,T,tcrit);
            flux_pr   = rainfall_1(P,T,tcrit);
            flux_qn   = melt_1(ddf,tcrit,T,S1,delta_t);
            flux_et1  = evap_7(S2,s2max,flux_epc,delta_t);
            flux_i    = interception_4(i_alpha,i_s,obj.t,tmax,flux_pr,delta_t);
            flux_q1f  = saturation_1(flux_pr+flux_qn,S2,s2max);
            flux_qw   = recharge_3(tw,S2);
            flux_et2  = evap_7(S3,se*s3max,flux_epc,delta_t);
            flux_q2f  = saturation_1(flux_qw,S3,s3max);
            flux_q2u  = baseflow_1(tu,S3);
            flux_qf   = baseflow_1(tc,S4);
            flux_qs   = baseflow_1(tc,S5);

            % stores ODEs
            dS1 = flux_ps - flux_qn;
            dS2 = flux_pr + flux_qn - flux_et1 - flux_i - flux_q1f - flux_qw;
            dS3 = flux_qw - flux_et2 - flux_q2f - flux_q2u;
            dS4 = flux_q1f + flux_q2f - flux_qf;
            dS5 = flux_q2u - flux_qs;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_epc, flux_ps,  flux_pr,  flux_qn,...
                      flux_et1, flux_i,   flux_q1f, flux_qw,...
                      flux_et2, flux_q2f, flux_q2u, flux_qf, flux_qs];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
        end
    end
end