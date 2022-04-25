classdef m_46_classic_12p_8s < MARRMoT_model
% Class for hydrologic conceptual model: CLASSIC

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Crooks, S. M., & Naden, P. S. (2007). CLASSIC: a semi-distributed 
% rainfall-runoff modelling system. Hydrology and Earth System Sciences, 
% 11(1), 516–531. http://doi.org/10.5194/hess-11-516-2007

    properties
        % model-specific attributes
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_46_classic_12p_8s()
            obj.numStores = 8;                                             % number of model stores
            obj.numFluxes = 21;                                            % number of model fluxes
            obj.numParams = 12;

            obj.JacobPattern  = [1,0,0,0,0,0,0,0;
                                 1,1,0,0,0,0,0,0;
                                 1,1,1,0,0,0,0,0;
                                 0,0,0,1,0,0,0,0;
                                 0,0,0,1,1,0,0,0;
                                 0,0,0,1,1,1,0,0;
                                 0,0,0,1,1,0,1,0;
                                 0,0,0,0,0,0,0,1];                         % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 1;       % fap, Fraction of catchment area that has permeable soils [-]
                                0.01, 0.99; % fdp, Fraction of depth of permeable soil that is store Px [-]
                                1, 2000;    % dp, Depth of permeable soil [mm]
                                0, 1;       % cq, Runoff coefficient for permeable soil [d-1]
                                0, 1;       % d1, Fraction of Ps that infiltrates into semi-permeable soil [-]
                                0, 1;       % fas, Fraction of (1-fap) that is fas [-]
                                0.01, 0.99; % fds, Fraction of depth of semi-permeable soil that is store Sx [-]
                                1, 2000;    % ds, Depth of semi-permeable soil [mm]
                                0, 1;       % d2, Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]
                                0, 1;       % cxq, Quick runoff coefficient for semi-permeable soil [d-1]
                                0, 1;       % cxs, Slow runoff coefficient for semi-permeable soil [d-1]
                                0, 1];      % cu, Runoff coefficient for impermeable soil [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5" "S6" "S7" "S8"};    % Names for the stores
            obj.FluxNames  = {"pp",  "ps",  "pi",  "epx", "ppx",...
                              "epy", "ppe", "q",   "psd", "psi",...
                              "esx", "psx", "esy", "pse", "psq",...
                              "pss", "xq",  "xs",  "ei",  "pie", "u"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = [4 6 11 13 19];                            % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 17 18 21];                              % Index or indices of fluxes to add to Streamflow
            obj.StoreSigns  = [1 -1 1 1 -1 1 1 1];                         % Signs to give to stores (-1 is a deficit store), only needed for water balance

        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta = obj.theta;
            fap   = theta(1);     % Fraction of catchment area that has permeable soils [-]
            fdp   = theta(2);     % Fraction of depth of permeable soil that is store Px [-]
            dp    = theta(3);     % Depth of permeable soil [mm]
            tf    = theta(6);     % Fraction of (1-fap) that is fas [-]
            fds   = theta(7);     % Fraction of depth of semi-permeable soil that is store Sx [-]
            ds    = theta(8);     % Depth of semi-permeable soil [mm]
            
            % auxiliary parameters;
            fas = (1-fap)*tf;     % Fraction of catchment area that has semi-permeable soils [-]
            fai = 1-fap-fas;      % Fraction of catchment area that has impermeable soils [-]
            pxm = fdp*dp;         % Depth of store Px [mm]
            pym = (1-fdp)*dp;     % Depth of store Py [mm]
            sxm = fds*ds;         % Depth of store Sx [mm]
            sym = (1-fds)*ds;     % Depth of store Sy [mm]
            obj.aux_theta = [fas, fai, pxm, pym, sxm, sym];
            
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            fap   = theta(1);     % Fraction of catchment area that has permeable soils [-]
            fdp   = theta(2);     % Fraction of depth of permeable soil that is store Px [-]
            dp    = theta(3);     % Depth of permeable soil [mm]
            cq    = theta(4);     % Runoff coefficient for permeable soil [d-1]
            d1    = theta(5);     % Fraction of Ps that infiltrates into semi-permeable soil [-]
            tf    = theta(6);     % Fraction of (1-fap) that is fas [-]
            fds   = theta(7);     % Fraction of depth of semi-permeable soil that is store Sx [-]
            ds    = theta(8);     % Depth of semi-permeable soil [mm]
            d2    = theta(9);     % Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]
            cxq   = theta(10);    % Quick runoff coefficient for semi-permeable soil [d-1]
            cxs   = theta(11);    % Slow runoff coefficient for semi-permeable soil [d-1]
            cu    = theta(12);    % Runoff coefficient for impermeable soil [d-1]
            
            % auxiliary parameters
            aux_theta = obj.aux_theta;
            fas = aux_theta(1);   % Fraction of catchment area that has semi-permeable soils [-]
            fai = aux_theta(2);   % Fraction of catchment area that has impermeable soils [-]
            pxm = aux_theta(3);   % Depth of store Px [mm]
            pym = aux_theta(4);   % Depth of store Py [mm]
            sxm = aux_theta(5);   % Depth of store Sx [mm]
            sym = aux_theta(6);   % Depth of store Sy [mm]
            
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
            S8 = S(8);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pp  = split_1(fap,P);
            flux_ps  = split_1(fas,P);
            flux_pi  = split_1(fai,P);
            flux_epx = evap_1(S1,fap*Ep,delta_t);
            flux_ppx = saturation_1(flux_pp,S1,pxm);
            flux_epy = evap_18(1.9,0.6523,pxm,S2+pxm,fap*Ep-flux_epx);
            flux_ppe = saturation_9(flux_ppx,S2,0.01);
            flux_q   = baseflow_1(cq,S3);
            flux_psd = split_1(1-d1,flux_ps);
            flux_psi = split_1(d1,flux_ps);
            flux_esx = evap_1(S4,fas*Ep,delta_t);
            flux_psx = saturation_1(flux_psi,S4,sxm);
            flux_esy = evap_18(1.9,0.6523,sxm,S4+sxm,fas*Ep-flux_esx);
            flux_pse = saturation_9(flux_psx,S5,0.01);
            flux_psq = split_1(d2,flux_pse+flux_psd);
            flux_pss = split_1(1-d2,flux_pse+flux_psd);
            flux_xq  = baseflow_1(cxq,S6);
            flux_xs  = baseflow_1(cxs,S7);
            flux_pie = effective_1(fai*flux_pi,0.5);
            flux_ei  = flux_pi - flux_pie;
            flux_u   = baseflow_1(cu,S8);

            % stores ODEs
            dS1 =    flux_pp  - flux_epx - flux_ppx;
            dS2 =  -(flux_ppx - flux_epy - flux_ppe);    
            dS3 =    flux_ppe - flux_q;    
            dS4 =    flux_psi - flux_esx - flux_psx;    
            dS5 =  -(flux_psx - flux_esy - flux_pse);    
            dS6 =    flux_psq - flux_xq;    
            dS7 =    flux_pss - flux_xs;    
            dS8 =    flux_pie - flux_u; 
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6 dS7 dS8];
            fluxes = [flux_pp,  flux_ps,  flux_pi,  flux_epx, flux_ppx,...
                      flux_epy, flux_ppe, flux_q,   flux_psd, flux_psi,...
                      flux_esx, flux_psx, flux_esy, flux_pse, flux_psq,...
                      flux_pss, flux_xq,  flux_xs,  flux_ei,  flux_pie, flux_u];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end