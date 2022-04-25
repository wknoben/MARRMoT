classdef m_05_ihacres_7p_1s < MARRMoT_model
% Class for hydrologic conceptual model: IHACRES

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Croke, B. F. W., & Jakeman, A. J. (2004). A catchment moisture deficit 
% module for the IHACRES rainfall-runoff model. Environmental Modelling and
% Software, 19(1), 1–5. http://doi.org/10.1016/j.envsoft.2003.09.001
%
% Littlewood, I. G., Down, K., Parker, J. R., & Post, D. A. (1997). IHACRES
% v1.0 User Guide.
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral 
% catchments. Water Resources Research, 33(1), 153–166. 
% http://doi.org/doi:10.1029/96WR02840

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_05_ihacres_7p_1s()
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 7;                                             % number of model fluxes
            obj.numParams = 7;
            
            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [1   , 2000;    % lp, Wilting point [mm]
                             1   , 2000;    % d, Threshold for flow generation [mm]
                             0   , 10;      % p, Flow response non-linearity [-]
                             0   , 1;       % alpha, Fast/slow flow division [-]
                             1   , 700;     % tau_q, Fast flow routing delay [d]
                             1   , 700;     % tau_s, Slow flow routing delay [d]
                             0   , 119];    % tau_d, flow delay [d]
            
            obj.StoreNames = {"S1"};                                       % Names for the stores
            obj.FluxNames  = {"Ea", "u", "uq", "us", "xq","xs", "Qt"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = 1;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = 7;                                         % Index or indices of fluxes to add to Streamflow
            obj.StoreSigns  = -1;                                          % Signs to give to stores (-1 is a deficit store), only needed for water balance
            
        end
        
        % INITialisation function
        function obj = init(obj)
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            tau_q   = theta(5);     % Fast flow routing delay [d]
            tau_s   = theta(6);     % Slow flow routing delay [d]
            tau_d   = theta(7);     % Pure time delay of total flow [d]
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh_q = uh_5_half(tau_q,delta_t);
            uh_s = uh_5_half(tau_s,delta_t);
            uh_t = uh_8_delay(tau_d,delta_t);
            
            obj.uhs = {uh_q, uh_s, uh_t};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            lp      = theta(1);     % Wilting point [mm]
            d       = theta(2);     % Threshold for flow generation [mm]
            p       = theta(3);     % Flow response non-linearity [-]
            alpha   = theta(4);     % Fast/slow flow division [-]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q = uhs{1};
            uh_s = uhs{2};
            uh_t = uhs{3};
            
            % stores
            S1 = S(1);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            
            % fluxes functions
            flux_ea   = evap_12(S1,lp,Ep);
            flux_u    = saturation_5(S1,d,p,P);
        	flux_uq   = split_1(alpha,flux_u);
            flux_us   = split_1(1-alpha,flux_u);            
            flux_xq   = route(flux_uq, uh_q);
            flux_xs   = route(flux_us, uh_s);
            flux_xt   = route(flux_xq + flux_xs, uh_t);

            % stores ODEs
            dS1 = -P + flux_ea + flux_u;
            
            % outputs
            dS = dS1;
            fluxes = [flux_ea, flux_u,  flux_uq, flux_us,...
                      flux_xq, flux_xs, flux_xt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q = uhs{1};
            uh_s = uhs{2};
            uh_t = uhs{3};
            
            % input fluxes to the unit hydrographs 
            fluxes = obj.fluxes(obj.t,:);
            flux_uq = fluxes(3);
            flux_us = fluxes(4);
            flux_xq = fluxes(5);
            flux_xs = fluxes(6);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            
            uh_q = update_uh(uh_q, flux_uq);
            uh_s = update_uh(uh_s, flux_us);
            uh_t = update_uh(uh_t, flux_xq + flux_xs);
            
            obj.uhs = {uh_q, uh_s, uh_t};
        end
    end
end