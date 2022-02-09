classdef m_14_topmodel_7p_2s < MARRMoT_model
% Class for hydrologic conceptual model: TOPMODEL

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Beven, K., Lamb, R., Quinn, P., Romanowicz, R., & Freer, J. (1995). 
% TOPMODEL. In V. P. Singh (Ed.), Computer Models of Watershed Hydrology 
% (pp. 627–668). Baton Rouge: Water Resources Publications, USA.
%
% Beven, K. J., & Kirkby, M. J. (1979). A physically based, variable 
% contributing area model of basin hydrology / Un modèle à base physique 
% de zone d’appel variable de l'hydrologie du bassin versant. Hydrological 
% Sciences Bulletin, 24(1), 43–69. http://doi.org/10.1080/02626667909491834
%
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. a., Vrugt, J. a., 
% Gupta, H. V., … Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between 
% hydrological models. Water Resources Research, 44(12). 
% http://doi.org/10.1029/2007WR006735
%
% Sivapalan, M., Beven, K., & Wood, E. F. (1987). On hydrologic similarity:
% 2. A scaled model of storm runoff production. Water Resources Research, 
% 23(12), 2266–2278. http://doi.org/10.1029/WR023i012p02266

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_14_topmodel_7p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 6;                                             % number of model fluxes
            obj.numParams = 7;
            
            obj.JacobPattern  = [1 1;...
                                 1 1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1, 2000;    % suzmax, Maximum soil moisture storage in unsatured zone [mm]
                             0.05, 0.95; % st, Threshold for flow generation and evap change as fraction of suzmax [-]
                             0,1;        % kd, Leakage to saturated zone flow coefficient [mm/d]
                             0.1,200;    % q0, Zero deficit base flow speed [mm/d]
                             0,1;        % m, Baseflow coefficient [mm-1]
                             1, 7.5;     % chi, Gamma distribution parameter [-]
                             0.1, 5];    % phi, Gamma distribution parameter [-]
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"qof", "peff", "ea", "qex", "qv", "qb"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = 3;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [1 4 6];                                   % Index or indices of fluxes to add to Streamflow
            obj.StoreSigns  = [1 -1];                                      % Signs to give to stores (-1 is a deficit store), only needed for water balance

        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            suzmax  = theta(1);     % Maximum soil moisture storage in unsatured zone [mm]
            st      = theta(2);     % Threshold for flow generation and evap change as fraction of suzmax [-]
            kd      = theta(3);     % Leakage to saturated zone flow coefficient [mm/d]
            q0      = theta(4);     % Zero deficit base flow speed [mm/d]
            f       = theta(5);     % Baseflow scaling coefficient [mm-1]
            chi     = theta(6);     % Gamma distribution parameter [-]
            phi     = theta(7);     % Gamma distribution parameter [-]
            mu      = 3;            % Gamma distribution parameter, fixed (Clark et al, 2008)
            lambda  = chi*phi+mu;   % Ac computation parameter, mean of the gamma distribution

            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_qof  = saturation_7(chi,phi,3,lambda,f,S2,P);
            flux_peff = P - flux_qof;
            flux_ea   = evap_3(st,S1,suzmax,Ep,delta_t);
            flux_qex  = saturation_1(flux_peff,S1,suzmax);
            flux_qv   = interflow_10(S1,kd,st*suzmax,suzmax-st*suzmax);
            flux_qb   = baseflow_4(q0,f,S2);

            % stores ODEs
            dS1 = flux_peff - flux_ea - flux_qex - flux_qv;
            dS2 = flux_qb - flux_qv;                                       % S2 is a deficit store
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_qof,  flux_peff, flux_ea,...
                      flux_qex, flux_qv, flux_qb];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end