classdef m_41_nam_10p_6s < MARRMoT_model
% Class for hydrologic conceptual model: NAM

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Nielsen, S. A., & Hansen, E. (1973). Numerical simulation of he rainfall-
% runoff process on a daily basis. Nordic Hydrology, (4), 171–190. 
% http://doi.org/https://doi.org/10.2166/nh.1973.0013

    properties
        % model-specific attributes
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_41_nam_10p_6s()
            obj.numStores = 6;                                             % number of model stores
            obj.numFluxes = 14;                                            % number of model fluxes
            obj.numParams = 10;
            
            obj.JacobPattern  = [1,0,0,0,0,0;
                                 1,1,1,0,0,0;
                                 1,1,1,0,0,0;
                                 1,1,1,1,0,0;
                                 0,1,1,0,1,0;
                                 1,1,1,0,0,1];                             % Jacobian matrix of model store ODEs
            
            obj.parRanges = [   0, 20;      % cs, Degree-day factor for snowmelt [mm/oC/d]
                                0, 1;       % cif, Runoff coefficient for interflow [d-1]
                                1, 2000;    % stot, Maximum total soil moisture depth [mm]
                                0, 0.99;    % cl1, Lower zone filling threshold for interflow generation [-]
                                0.01, 0.99; % f1, Fraction of total soil depth that makes up lstar
                                0, 1;       % cof, Runoff coefficient for overland flow [d-1]
                                0, 0.99;    % cl2, Lower zone filling threshold for overland flow generation [-]
                                0, 1;       % k0, Overland flow routing delay [d-1]
                                0, 1;       % k1, Interflow routing delay [d-1]
                                0, 1];      % kb, Baseflow routing delay [d-1]
                            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5" "S6"};              % Names for the stores
            obj.FluxNames  = {"ps", "pr", "m",  "eu", "pn", "of", "inf",...
                              "if", "dl", "gw", "el", "qo", "qi", "qb"};   % Names for the fluxes
            
            obj.FluxGroups.Ea = [4 11];                                    % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [12 13 14];                                % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INIT is run automatically as soon as both theta and delta_t are
        % set (it is therefore ran only once at the beginning of the run. 
        % Use it to initialise all the model parameters (in case there are
        % derived parameters) and unit hydrographs and set minima and
        % maxima for stores based on parameters.
        function obj = init(obj)
            % parameteres
            theta = obj.theta;
            stot  = theta(3);     % Maximum total soil moisture depth [mm]
            fl    = theta(5);     % Fraction of total soil depth that makes up lstar
            
            % set auxiliary parameters
            lstar = fl*stot;      % Maximum lower zone storage [mm]
            ustar = (1-fl)*stot;  % Upper zone maximum storage [mm]
            obj.aux_theta = [lstar; ustar];
            
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            cs    = theta(1);     % Degree-day factor for snowmelt [mm/oC/d]
            cif   = theta(2);     % Runoff coefficient for interflow [d-1]
            cl1   = theta(4);     % Lower zone filling threshold for interflow generation [-]
            cof   = theta(6);     % Runoff coefficient for overland flow [d-1]
            cl2   = theta(7);     % Lower zone filling threshold for overland flow generation [-]
            k0    = theta(8);     % Overland flow routing delay [d-1]
            k1    = theta(9);     % Interflow routing delay [d-1]
            kb    = theta(10);    % Baseflow routing delay [d-1]
            
            % auxiliary parameters
            aux_theta = obj.aux_theta;
            lstar = aux_theta(1);  % Maximum lower zone storage [mm]
            ustar = aux_theta(2);  % Upper zone maximum storage [mm]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            S6 = S(6);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ps   = snowfall_1(P,T,0);
            flux_pr   = rainfall_1(P,T,0);
            flux_m    = melt_1(cs,0,T,S1,delta_t);
            flux_eu   = evap_1(S2,Ep,delta_t);
            flux_pn   = saturation_1(flux_pr+flux_m,S2,ustar);
            flux_of   = interflow_6(cof,cl2,flux_pn,S3,lstar);
            flux_inf  = flux_pn-flux_of;
            flux_if   = interflow_6(cif,cl1,S2,S3,lstar);
            flux_dl   = split_1(1-S3/lstar,flux_inf);
            flux_gw   = split_1(S3/lstar,flux_inf);
            flux_el   = evap_15(Ep,S3,lstar,S2,0.01,delta_t);
            flux_qo   = interflow_5(k0,S4);
            flux_qi   = interflow_5(k1,S5);
            flux_qb   = baseflow_1(kb,S6);
            
            % stores ODEs
            dS1 = flux_ps - flux_m;
            dS2 = flux_pr + flux_m - flux_eu - flux_if - flux_pn;    
            dS3 = flux_dl - flux_el;    
            dS4 = flux_of - flux_qo;    
            dS5 = flux_if - flux_qi;    
            dS6 = flux_gw - flux_qb; 
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6];
            fluxes = [flux_ps, flux_pr,  flux_m,  flux_eu, flux_pn,...
                      flux_of, flux_inf, flux_if, flux_dl, flux_gw,...
                      flux_el, flux_qo,  flux_qi, flux_qb];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
        end
    end
end