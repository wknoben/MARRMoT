classdef m_22_vic_10p_3s < MARRMoT_model
% Class for hydrologic conceptual model:  VIC

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. A., Vrugt, J. A.,
% Gupta, H. V., � Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between
% hydrological models. Water Resources Research, 44(12), W00B02.
% http://doi.org/10.1029/2007WR006735
%
% Liang, X., Lettenmaier, D. P., Wood, E. F., & Burges, S. J. (1994). A
% simple hydrologically based model of land surface water and energy fluxes
% for general circulation models. Journal of Geophysical Research, 99,
% 14415�14428.

    properties
        % model-specific attributes
        aux_theta     % auxiliary parameter set
    end
    methods

        % creator method
        function obj = m_22_vic_10p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 11;                                            % number of model fluxes
            obj.numParams = 10;

            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs

            obj.parRanges = [   0.1 , 5;        % ibar,     Mean interception capacity [mm]
                                0   , 1;        % idelta,   Seasonal interception change as fraction of mean [-]
                                1   , 365;      % ishift,   Maximum interception peak timing [-]
                                1   , 2000;     % stot,    Maximum soil moisture capacity [mm]
                                0.01, 0.99;     % fsm, Fraction of stot that constitutes maximum soil moisture smmax [-]
                                0   ,10;        % b,        Infiltration excess shape parameter [-]
                                0   , 1;        % k1,       Percolation time parameter [d-1]
                                0   ,10;        % c1,       Percolation non-linearity parameter [-]
                                0   ,1;         % k2,       Baseflow time parameter [d-1]
                                1   , 5];       % c2,       Baseflow non-linearity parameter

            obj.StoreNames = {"S1", "S2" "S3"};                             % Names for the stores
            obj.FluxNames  = {"ei",   "peff", "iex", "qie",  "inf", "et1",...
                              "qex1", "pc",   "et2", "qex2", "qb"};        % Names for the fluxes

            obj.FluxGroups.Ea = [1 6 9];                                   % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [4 7 10 11];                               % Index or indices of fluxes to add to Streamflow

        end

        % INITialisation function
        function obj = init(obj)
            %Parameters
            theta = obj.theta;
            stot    = theta(4);     % Total available storage [mm]
            fsm     = theta(5);     % Fraction of stot that constitutes maximum soil mositure storage [-]

            % Auxiliary parameter
            smmax   = fsm*stot;     % Maximum soil moisture capacity [mm]
            gwmax   = (1-fsm)*stot; % Maximum groundwater storage [mm]
            tmax    = 365.25;       % Length of one growing cycle [d]
            obj.aux_theta = [smmax; gwmax; tmax];

        end

        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            ibar    = theta(1);     % Mean interception capacity [mm]
            idelta  = theta(2);     % Seasonal interception change as fraction of mean [-]
            ishift  = theta(3);     % Maximum interception peak timing [-]
            b       = theta(6);     % Infiltration excess shape parameter [-]
            k1      = theta(7);     % Percolation time parameter [d-1]
            c1      = theta(8);     % Percolation non-linearity parameter [-]
            k2      = theta(9);     % Baseflow time parameter [d-1]
            c2      = theta(10);    % Baseflow non-linearity parameter

            aux_theta = obj.aux_theta;
            smmax = aux_theta(1);
            gwmax = aux_theta(2);
            tmax = aux_theta(3);

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
            T  = climate_in(3);

            % fluxes functions
            aux_imax  = phenology_2(ibar,idelta,ishift,obj.t,tmax,delta_t);
            flux_ei   = evap_7(S1,aux_imax,Ep,delta_t);
            flux_peff = interception_1(P,S1,aux_imax);
            flux_iex  = excess_1(S1,aux_imax,delta_t);
            flux_qie  = saturation_2(S2,smmax,b,flux_peff+flux_iex);
            flux_inf  = effective_1(flux_peff+flux_iex,flux_qie);
            flux_et1  = evap_7(S2,smmax,max(0,Ep-flux_ei),delta_t);
            flux_qex1 = saturation_1(flux_inf,S2,smmax);
            flux_pc   = percolation_5(k1,c1,S2,smmax,delta_t);
            flux_et2  = evap_7(S3,gwmax,max(0,Ep-flux_ei-flux_et1),delta_t);
            flux_qex2 = saturation_1(flux_pc,S3,gwmax);
            flux_qb   = baseflow_5(k2,c2,S3,gwmax,delta_t);

            % stores ODEs
            dS1 = P - flux_ei - flux_peff - flux_iex;
            dS2 = flux_inf - flux_et1 - flux_qex1 - flux_pc;
            dS3 = flux_pc  - flux_et2 - flux_qex2 - flux_qb;

            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ei  flux_peff flux_iex  flux_qie ...
                      flux_inf flux_et1  flux_qex1 flux_pc ...
                      flux_et2 flux_qex2 flux_qb];
        end

        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end
