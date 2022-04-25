classdef m_25_tcm_6p_4s < MARRMoT_model
% Class for hydrologic conceptual model: Thames Catchment Model

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Moore, R. J., & Bell, V. A. (2001). Comparison of rainfall-runoff models
% for flood forecasting. Part 1: Literature review of models. Bristol:
% Environment Agency.

    properties
        % model-specific attributes
        aux_theta     % auxiliary parameter set
    end
    methods

        % creator method
        function obj = m_25_tcm_6p_4s()
            obj.numStores = 4;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 6;

            obj.JacobPattern  = [1,0,0,0;
                                 1,1,0,0;
                                 1,1,1,0
                                 0,0,1,1];                                 % Jacobian matrix of model store ODEs

            obj.parRanges = [   0, 1;     % phi, Fraction preferential recharge [-]
                                1, 2000;  % rc, Maximum soil moisture depth [mm]
                                0, 1;     % gam, Fraction of Ep reduction with depth [-]
                                0, 1;     % k1, Runoff coefficient [d-1]
                                0, 1;     % fa, Fraction of mean(P) that forms abstraction rate [mm/d]
                                0, 1];    % k2, Runoff coefficient [mm-1 d-1]

            obj.StoreNames = {"S1", "S2" "S3" "S4"};                        % Names for the stores
            obj.FluxNames  = {"pn", "en",   "pby",   "pin",  "ea",...
                              "et", "qex1", "qex2", "quz", "a", "q"};      % Names for the fluxes

            obj.FluxGroups.Ea = [2 5 6];                                   % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 11;                                        % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.Abstraction = 10;                               % Index or abstraction flux (just needed for water balance)
            obj.StoreSigns  = [1 -1 1 1];                                  % Signs to give to stores (-1 is a deficit store), only needed for water balance
        end

        % INITialisation function
        function obj = init(obj)
            fa = obj.theta(5);    % Fraction of average P abstracted [-]
            P  = obj.input_climate(:,1);

            ca = fa * mean(P);    % Abstraction rate [mm/day]
            obj.aux_theta(1) = ca;
        end

        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            phi   = theta(1);     % Fraction preferential recharge [-]
            rc    = theta(2);     % Maximum soil moisture depth [mm]
            gam   = theta(3);     % Fraction of Ep reduction with depth [-]
            k1    = theta(4);     % Runoff coefficient [d-1]
            k2    = theta(6);     % Runoff coefficient [mm-1 d-1]

            ca    = obj.aux_theta(1);    % Abstraction rate [mm/day]

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
            flux_pn     = effective_1(P,Ep);
            flux_en     = P - flux_pn;
            flux_pby    = split_1(phi,flux_pn);
            flux_pin    = split_1(1-phi,flux_pn);
            flux_ea     = evap_1(S1,Ep,delta_t);
            flux_et     = evap_16(gam,Inf,S1,0.01,Ep,delta_t);
            flux_qex1   = saturation_1(flux_pin,S1,rc);
            flux_qex2   = saturation_9(flux_qex1,S2,0.01);
            flux_quz    = baseflow_1(k1,S3);
            flux_a      = abstraction_1(ca);
            flux_q      = baseflow_6(k2,0,S4,delta_t);

            % stores ODEs
            dS1 = flux_pin  - flux_ea   - flux_qex1;
            dS2 = flux_et   + flux_qex2 - flux_qex1;
            dS3 = flux_qex2 + flux_pby  - flux_quz;
            dS4 = flux_quz  - flux_a    - flux_q;

            % outputs
            dS = [dS1 dS2 dS3 dS4];
            fluxes = [flux_pn, flux_en, flux_pby,  flux_pin,...
                      flux_ea, flux_et, flux_qex1, flux_qex2,...
                      flux_quz, flux_a, flux_q];
        end

        % STEP runs at the end of every timestep.
        function obj = step(obj)
        end
    end
end
