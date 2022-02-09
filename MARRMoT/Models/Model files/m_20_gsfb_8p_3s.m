classdef m_20_gsfb_8p_3s < MARRMoT_model
% Class for hydrologic conceptual model: GSFB

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Nathan, R. J., & McMahon, T. A. (1990). SFB model part l . Validation of 
% fixed model parameters. In Civil Eng. Trans. (pp. 157–161).
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral catchments.
% Water Resources Research, 33(1), 153–166. http://doi.org/doi:10.1029/96WR02840

    properties
        % model-specific attributes
    end
    methods
        
        % creator method

        function obj = m_20_gsfb_8p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 6;                                             % number of model fluxes
            obj.numParams = 8; 

            obj.JacobPattern  = [1,0,1;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 1;           % c, Recharge time coeffcient [d-1]
                                0.05, 0.95;     % ndc, Threshold fraction of Smax [-]
                                1, 2000;        % smax, Maximum soil moisture storage [mm]
                                0, 20;          % emax, Maximum evaporation flux [mm/d]
                                0, 200;         % frate, Maximum infiltration rate [mm/d]
                                0, 1;           % b, Fraction of subsurface flow that is baseflow [-]
                                0, 1;           % dpf, Baseflow time coefficient [d-1]
                                1, 300];        % sdrmax, Threshold before baseflow can occur [mm]
            
            obj.StoreNames = {"S1", "S2" "S3"};                             % Names for the stores
            obj.FluxNames  = {"ea", "qs", "f", "qb", "dp", "qdr"};         % Names for the fluxes
            
            obj.FluxGroups.Ea = 1;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [2 4];                                     % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            c       = theta(1);     % Recharge time coeffcient [d-1]
            ndc     = theta(2);     % Threshold fraction of Smax [-]
            smax    = theta(3);     % Maximum soil moisture storage [mm]
            emax    = theta(4);     % Maximum evaporation flux [mm/d]
            frate   = theta(5);     % Maximum infiltration rate [mm/d]
            b       = theta(6);     % Fraction of subsurface flow that is baseflow [-]
            dpf     = theta(7);     % Baseflow time coefficient [d-1]
            sdrmax  = theta(8);     % Threshold before baseflow can occur [mm]
            
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
            flux_ea  = evap_20(emax,ndc,S1,smax,Ep,delta_t);
            flux_qs  = saturation_1(P,S1,smax);
            flux_f   = interflow_11(frate,ndc*smax,S1,delta_t);
            flux_qb  = baseflow_9(b*dpf,sdrmax,S2);
            flux_dp  = baseflow_1((1-b)*dpf,S2);
            flux_qdr = recharge_5(c,ndc*smax,S3,S1);

            % stores ODEs
            dS1 = P + flux_qdr - flux_ea - flux_qs - flux_f;
            dS2 = flux_f - flux_qb - flux_dp;    
            dS3 = flux_dp - flux_qdr;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ea, flux_qs, flux_f, flux_qb, flux_dp, flux_qdr];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end