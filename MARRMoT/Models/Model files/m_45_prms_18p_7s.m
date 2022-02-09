classdef m_45_prms_18p_7s < MARRMoT_model
% Class for hydrologic conceptual model: PRMS

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Leavesley, G. H., R. Lichty, B. Troutman, and L. Saindon (1983), 
% Precipitation-Runo Modeling System: User's Manual. U.S. Geological 
% Survey, Water-Resources Investigations Report 83-4238, 207
%
% Markstrom, S. L., S. Regan, L. E. Hay, R. J. Viger, R. M. T. Webb, R. A. 
% Payn, and J. H. LaFontaine (2015), PRMS-IV, the Precipitation-Runoff 
% Modeling System, Version 4. In U.S. Geological Survey Techniques and
% Methods, book 6, chap. B7, 158. doi: http://dx.doi.org/10.3133/tm6B7

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_45_prms_18p_7s()
            obj.numStores = 7;                                             % number of model stores
            obj.numFluxes = 25;                                            % number of model fluxes
            obj.numParams = 18; 

            obj.JacobPattern  = [1,0,0,0,0,0,0;
                                 0,1,0,0,0,0,0;
                                 1,0,1,0,0,0,0;
                                 1,1,1,1,0,0,0;
                                 1,1,1,1,1,0,0;
                                 1,1,1,1,1,1,0;
                                 1,1,1,1,1,1,1];                           % Jacobian matrix of model store ODEs
            
            obj.parRanges = [-3,  5;         % tt, Temperature threshold for snowfall and melt [oC]
                              0, 20;         % ddf,  Degree-day factor for snowmelt [mm/oC/d]
                              0,  1;         % alpha, Fraction of rainfall on soil moisture going to interception [-] 
                              0,  1;         % beta, Fraction of catchment where rain goes to soil moisture [-]
                              0,  5;         % stor, Maximum interception capcity [mm]
                              0, 50;         % retip, Maximum impervious area storage [mm]
                              0,  1;         % fscn, Fraction of SCX where SCN is located [-]
                              0,  1;         % scx, Maximum contributing fraction area to saturation excess flow [-]
                              0.005, 0.995;  % flz, Fraction of total soil moisture that is the lower zone [-]
                              1, 2000;       % stot, Total soil moisture storage [mm]: REMX+SMAX
                              0, 20;         % cgw, Constant drainage to deep groundwater [mm/d]
                              1, 300;        % resmax, Maximum flow routing reservoir storage (used for scaling only, there is no overflow) [mm]
                              0,  1;         % k1, Groundwater drainage coefficient [d-1]
                              1,  5;         % k2, Groundwater drainage non-linearity [-]
                              0,  1;         % k3, Interflow coefficient 1 [d-1]
                              0,  1;         % k4, Interflow coefficient 2 [mm-1 d-1]
                              0,  1;         % k5, Baseflow coefficient [d-1]
                              0,  1];        % k6, Groundwater sink coefficient [d-1]

            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5" "S6" "S7"};         % Names for the stores
            obj.FluxNames  = {"ps",   "pr",  "pim", "psm", "pby"...
                              "pin",  "ptf", "m",   "mim", "msm"...
                              "sas",  "sro", "inf", "pc",  "excs"...
                              "qres", "sep", "gad", "ras", "bas"...
                              "snk",  "ein", "eim", "ea",  "et"};          % Names for the fluxes
            
            obj.FluxGroups.Ea = [22 23 24 25];                             % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [11 12 19 20];                             % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.Sink  = 21;                                     % index of sink fluxes

        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            tt      = theta(1);     % Temperature threshold for snowfall and melt [oC]
            ddf     = theta(2);     % Degree-day factor for snowmelt [mm/oC/d]
            alpha   = theta(3);     % Fraction of rainfall on soil moisture going to interception [-] 
            beta    = theta(4);     % Fraction of catchment where rain goes to soil moisture [-]
            stor    = theta(5);     % Maximum interception capcity [mm]
            retip   = theta(6);     % Maximum impervious area storage [mm]
            fscn    = theta(7);     % Fraction of SCX where SCN is located [-]
            scx     = theta(8);     % Maximum contributing fraction area to saturation excess flow [-]
            scn     = fscn*scx;     % Minimum contributing fraction area to saturation excess flow [-]
            flz     = theta(9);     % Fraction of total soil moisture that is the lower zone [-]
            stot    = theta(10);    % Total soil moisture storage [mm]: REMX+SMAX
            remx    = (1-flz)*stot; % Maximum upper soil moisture storage [mm]
            smax    = flz*stot;     % Maximum lower soil moisture storage [mm] 
            cgw     = theta(11);    % Constant drainage to deep groundwater [mm/d]
            resmax  = theta(12);    % Maximum flow routing reservoir storage (used for scaling only, there is no overflow) [mm]
            k1      = theta(13);    % Groundwater drainage coefficient [d-1]
            k2      = theta(14);    % Groundwater drainage non-linearity [-]
            k3      = theta(15);    % Interflow coefficient 1 [d-1]
            k4      = theta(16);    % Interflow coefficient 2 [mm-1 d-1]
            k5      = theta(17);    % Baseflow coefficient [d-1]
            k6      = theta(18);    % Groundwater sink coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            S6 = S(6);
            S7 = S(7);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ps  = snowfall_1(P,T,tt);
            flux_pr  = rainfall_1(P,T,tt);
            flux_pim = split_1(1-beta,flux_pr);
            flux_psm = split_1(beta,flux_pr);
            flux_pby = split_1(1-alpha,flux_psm);
            flux_pin = split_1(alpha,flux_psm);
            flux_ptf = interception_1(flux_pin,S2,stor);
            flux_m   = melt_1(ddf,tt,T,S1,delta_t);
            flux_mim = split_1(1-beta,flux_m);
            flux_msm = split_1(beta,flux_m);
            flux_sas = saturation_1(flux_pim+flux_mim,S3,retip);
            flux_sro = saturation_8(scn,scx,S4,remx,flux_msm+flux_ptf+flux_pby);
            flux_inf = effective_1(flux_msm+flux_ptf+flux_pby,flux_sro);
            flux_pc  = saturation_1(flux_inf,S4,remx);
            flux_excs= saturation_1(flux_pc,S5,smax);
            flux_sep = recharge_7(cgw,flux_excs);
            flux_qres= effective_1(flux_excs,flux_sep);
            flux_gad = recharge_2(k2,S6,resmax,k1);
            flux_ras = interflow_4(k3,k4,S6);
            flux_bas = baseflow_1(k5,S7);
            flux_snk = baseflow_1(k6,S7);
            flux_ein = evap_1(S2,beta*Ep,delta_t);
            flux_eim = evap_1(S3,(1-beta)*Ep,delta_t);
            flux_ea  = evap_7(S4,remx,Ep-flux_ein-flux_eim,delta_t);
            flux_et  = evap_15(Ep-flux_ein-flux_eim-flux_ea,S5,smax,S4,Ep-flux_ein-flux_eim,delta_t);

            % stores ODEs
            dS1 = flux_ps  - flux_m;
            dS2 = flux_pin - flux_ein - flux_ptf;    
            dS3 = flux_pim + flux_mim - flux_eim - flux_sas;
            dS4 = flux_inf - flux_ea  - flux_pc;
            dS5 = flux_pc  - flux_et  - flux_excs;
            dS6 = flux_qres- flux_gad - flux_ras;
            dS7 = flux_sep + flux_gad - flux_bas - flux_snk;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6 dS7];
            fluxes = [flux_ps,   flux_pr,  flux_pim, flux_psm, flux_pby...
                      flux_pin,  flux_ptf, flux_m,   flux_mim, flux_msm...
                      flux_sas,  flux_sro, flux_inf, flux_pc,  flux_excs...
                      flux_qres, flux_sep, flux_gad, flux_ras, flux_bas...
                      flux_snk,  flux_ein, flux_eim, flux_ea,  flux_et];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
        end
    end
end