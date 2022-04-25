classdef m_42_hycymodel_12p_6s < MARRMoT_model
% Class for hydrologic conceptual model: HYCYMODEL

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Fukushima, Y. (1988). A model of river flow forecasting for a small 
% forested mountain catchment. Hydrological Processes, 2(2), 167–185.

    properties
        % model-specific attributes
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_42_hycymodel_12p_6s()
            obj.numStores = 6;                                             % number of model stores
            obj.numFluxes = 18;                                            % number of model fluxes
            obj.numParams = 12; 

            obj.JacobPattern  = [1,0,0,0,0,0;
                                 1,1,0,0,0,0;
                                 1,1,1,0,0,0;
                                 0,0,1,1,0,0;
                                 1,1,1,0,1,0;
                                 0,0,0,0,0,1];                             % Jacobian matrix of model store ODEs
            
            obj.parRanges = [   0, 1;           % c,    Fraction area that is channel [-]
                                0, 5;           % imax, Maximum total interception storage [mm]
                                0, 1;           % a,    Fraction stem/trunk interception [-]
                                0.01, 0.99;     % fi2,  Fraction of total interception that is trunk/stem interception [mm]
                                0, 1;           % kin,  Infiltration runoff coefficient [d-1]
                                1, 2000;        % D50,  Soil depth where 50% of area contributes to effective flow [mm]
                                0.01, 0.99;     % fd16, Soil depth where 16% of area contributes to effective flow [mm]
                                1, 2000;        % sbc, Soil depth where evaporation rate starts to decline [mm]
                                0, 1;           % kb, Baseflow runoff coefficient [d-1]
                                1, 5;           % pb, Baseflow non-linearity [-]
                                0, 1;           % kh, Hillslope runoff coefficient [d-1]
                                0, 1];          % kc, Channel runoff coefficient [d-1]
                             
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5" "S6"};              % Names for the stores
            obj.FluxNames  = {"rc",  "rg", "eic", "qie", "qis", "rt"...
                              "eis", "rs", "rn",  "esu", "re",  "qin"...
                              "esb", "qb", "qh",  "qc",  "ec",  "qt"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = [3 7 10 13 17];                            % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [18];                                      % Index or indices of fluxes to add to Streamflow

        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta = obj.theta;
            imax    = theta(2);     % Maximum total interception storage [mm]
            fi2     = theta(4);     % Fraction of total interception that is trunk/stem interception [mm]
            d50     = theta(6);     % Soil depth where 50% of area contributes to effective flow [mm]
            fd16    = theta(7);     % Fraction of D50 that is D16 [-]

            % Auxiliary parameters 
            i1max   = (1-fi2)*imax; % Maximum canopy interception [mm]
            i2max   = fi2*imax;     % Maximum trunk/stem interception [mm]
            d16     = fd16*d50;     % Soil depth where 16% of area contributes to effective flow [mm]
            obj.aux_theta = [i1max, i2max, d16];

        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            c       = theta(1);     % Fraction area that is channel [-]
            imax    = theta(2);     % Maximum total interception storage [mm]
            a       = theta(3);     % Fraction stem/trunk interception [-]
            fi2     = theta(4);     % Fraction of total interception that is trunk/stem interception [mm]
            kin     = theta(5);     % Infiltration runoff coefficient [d-1]
            d50     = theta(6);     % Soil depth where 50% of area contributes to effective flow [mm]
            fd16    = theta(7);     % Fraction of D50 that is D16 [-]
            sbc     = theta(8);     % Soil depth where evaporation rate starts to decline [mm]
            kb      = theta(9);     % Baseflow runoff coefficient [d-1]
            pb      = theta(10);    % Baseflow non-linearity [-]
            kh      = theta(11);    % Hillslope runoff coefficient [d-1]
            kc      = theta(12);    % Channel runoff coefficient [d-1]

            % auxiliary parameters
            aux_theta = obj.aux_theta;
            i1max   = aux_theta(1); % Maximum canopy interception [mm]
            i2max   = aux_theta(2); % Maximum trunk/stem interception [mm]
            d16     = aux_theta(3); % Soil depth where 16% of area contributes to effective flow [mm]
            
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
            flux_rc  = split_1(c,P);
            flux_rg  = split_1(1-c,P);
            flux_eic = evap_1(S1,(1-c)*Ep,delta_t);
            flux_qie = interception_1(flux_rg,S1,i1max);
            flux_qis = split_1(a,flux_qie);
            flux_rt  = split_1(1-a,flux_qie);
            flux_eis = evap_1(S2,(1-c)*Ep,delta_t);
            flux_rs  = interception_1(flux_qis,S2,i2max);
            flux_rn  = flux_rt + flux_rs;
            flux_esu = evap_1(S3,(1-c)*Ep,delta_t);
            flux_re  = saturation_13(d50,d16,S3,flux_rn);
            flux_qin = recharge_3(kin,S3);
            flux_esb = evap_3(1,S4,sbc,max(0,(1-c)*Ep-flux_esu),delta_t);
            flux_qb  = baseflow_7(kb,pb,S4,delta_t);
            flux_qh  = interflow_3(kh,5/3,S5,delta_t);
            flux_qc  = interflow_3(kc,5/3,S6,delta_t);
            flux_qt  = effective_1(flux_qb+flux_qh+flux_qc,c*Ep);
            flux_ec  = (flux_qb+flux_qh+flux_qc) - flux_qt;

            % stores ODEs
            dS1 = flux_rg  - flux_eic - flux_qie;
            dS2 = flux_qis - flux_eis - flux_rs;
            dS3 = flux_rn  - flux_re  - flux_esu - flux_qin;
            dS4 = flux_qin - flux_esb - flux_qb;
            dS5 = flux_re  - flux_qh;
            dS6 = flux_rc  - flux_qc;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6];
            fluxes = [flux_rc,  flux_rg, flux_eic, flux_qie, flux_qis, flux_rt...
                      flux_eis, flux_rs, flux_rn,  flux_esu, flux_re,  flux_qin...
                      flux_esb, flux_qb, flux_qh,  flux_qc,  flux_ec,  flux_qt];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
        end
    end
end