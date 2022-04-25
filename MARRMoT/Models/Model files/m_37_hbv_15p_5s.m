classdef m_37_hbv_15p_5s < MARRMoT_model
% Class for hydrologic conceptual model: HBV-96

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Lindström, G., Johansson, B., Persson, M., Gardelin, M., & Bergström, S. 
% (1997). Development and test of the distributed HBV-96 hydrological model. 
% Journal of Hydrology, 201, 272–288. 
% https://doi.org/https://doi.org/10.1016/S0022-1694(97)00041-3

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_37_hbv_15p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 13;                                            % number of model fluxes
            obj.numParams = 15;
            
            obj.JacobPattern  = [1,1,0,0,0;
                                 1,1,0,0,0;
                                 1,1,1,1,0;
                                 1,1,1,1,0;
                                 0,0,0,1,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [-3 ,   5;      % TT, threshold temperature for snowfall [oC] 
                              0 ,   17;     % TTI, interval length of rain-snow spectrum [oC]
                             -3 ,   3;      % TTM, threshold temperature for snowmelt [oC]
                              0 ,   1;      % CFR, coefficient of refreezing of melted snow [-]
                              0,   20;      % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
                              0 ,   1;      % WHC, maximum water holding content of snow pack [-]
                              0 ,   4;      % CFLUX, maximum rate of capillary rise [mm/d]
                              1 ,   2000;   % FC, maximum soil moisture storage [mm]
                              0.05,    0.95;% LP, wilting point as fraction of FC [-]
                              0 ,   10;     % BETA, non-linearity coefficient of upper zone recharge [-]
                              0 ,   1;      % K0, runoff coefficient from upper zone [d-1] 
                              0 ,   4;      % ALPHA, non-linearity coefficient of runoff from upper zone [-]
                              0 ,   20;     % PERC, maximum rate of percolation to lower zone [mm/d]
                              0 ,   1;      % K1, runoff coefficient from lower zone [d-1]
                              1 ,   120];   % MAXBAS, flow routing delay [d]                 
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"sf", "refr", "melt", "rf",   "in", "se", "cf",...
                              "ea", "r",    "q0",   "perc", "q1", "qt"};   % Names for the fluxes
            
            obj.FluxGroups.Ea = 8;                                           % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 13;                                          % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            %parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            maxbas  = theta(15);    % MAXBAS, flow routing delay [d]
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(maxbas,delta_t);
            obj.uhs = {uh};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            tt      = theta(1);     % TT, middle of snow-rain interval [oC]
            tti     = theta(2);     % TTI, interval length of rain-snow spectrum [oC]
            ttm     = theta(3);     % TTM, threshold temperature for snowmelt [oC]
            cfr     = theta(4);     % CFR, coefficient of refreezing of melted snow [-]
            cfmax   = theta(5);     % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
            whc     = theta(6);     % WHC, maximum water holding content of snow pack [-]
            cflux   = theta(7);   	% CFLUX, maximum rate of capillary rise [mm/d]
            fc      = theta(8);     % FC, maximum soil moisture storage [mm]
            lp      = theta(9);     % LP, wilting point as fraction of FC [-]
            beta    = theta(10);    % BETA, non-linearity coefficient of upper zone recharge [-]
            k0      = theta(11);    % K0, runoff coefficient from upper zone [d-1], 
            alpha   = theta(12);    % ALPHA, non-linearity coefficient of runoff from upper zone [-]
            perc    = theta(13);    % PERC, maximum rate of percolation to lower zone [mm/d]
            k1      = theta(14);    % K1, runoff coefficient from lower zone [d-1]
            maxbas  = theta(15);    % MAXBAS, flow routing delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % stores at previous timestep
            t = obj.t;                             % this time step
            if t == 1
                S2old=obj.S0(2);
            else
                S2old = obj.stores(t-1,2);
            end
            
            % climate input
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_sf   = snowfall_2(P,T,tt,tti);
            flux_refr = refreeze_1(cfr,cfmax,ttm,T,S2,delta_t);
            flux_melt = melt_1(cfmax,ttm,T,S1,delta_t);
            flux_rf   = rainfall_2(P,T,tt,tti);
            flux_in   = infiltration_3(flux_rf+flux_melt,S2,whc*S1);
            flux_se   = excess_1(S2old,whc*S1,delta_t);
            flux_cf   = capillary_1(cflux,S3,fc,S4,delta_t);
            flux_ea   = evap_3(lp,S3,fc,Ep,delta_t);
            flux_r    = recharge_2(beta,S3,fc,flux_in+flux_se);
            flux_q0   = interflow_2(k0,S4,alpha,delta_t);
            flux_perc = percolation_1(perc,S4,delta_t);
            flux_q1   = baseflow_1(k1,S5);
            flux_qt   = route(flux_q0 + flux_q1, uh);

            % stores ODEs
            dS1 = flux_sf   + flux_refr - flux_melt;
            dS2 = flux_rf   + flux_melt - flux_refr - flux_in - flux_se;
            dS3 = flux_in   + flux_se   + flux_cf   - flux_ea - flux_r;
            dS4 = flux_r    - flux_cf   - flux_q0   - flux_perc;
            dS5 = flux_perc - flux_q1;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_sf,   flux_refr, flux_melt, flux_rf, flux_in,...
                      flux_se,   flux_cf,   flux_ea,   flux_r,  flux_q0,...
                      flux_perc, flux_q1, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % input fluxes to the unit hydrographs
            fluxes = obj.fluxes(obj.t,:);
            flux_q0 = fluxes(10);
            flux_q1 = fluxes(12);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh = update_uh(uh, flux_q0 + flux_q1);
            obj.uhs = {uh};
        end
    end
end