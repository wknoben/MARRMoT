classdef m_07_gr4j_4p_2s < MARRMoT_model
% Class for hydrologic conceptual model: GR4J

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Perrin, C., Michel, C., & Andréassian, V. (2003). Improvement of a 
% parsimonious model for streamflow simulation. Journal of Hydrology, 
% 279(1-4), 275–289. http://doi.org/10.1016/S0022-1694(03)00225-7
%
% Santos, L., Thirel, G., & Perrin, C. (2017). State-space representation 
% of a bucket-type rainfall-runoff model: a case study with State-Space GR4
% (version 1.0). Geoscientific Model Development Discussions, 1–22. 
% http://doi.org/10.5194/gmd-2017-264

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_07_gr4j_4p_2s()            
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 13;                                            % number of model fluxes
            obj.numParams = 4;
            
            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1   , 2000;    % x1 [mm]
                            -20  , 20;      % x2 [mm/d]
                             1   , 300;     % x3 [mm]
                             .5  , 15];     % x4 [d]
            
            obj.StoreNames = {"S1", "S2"};                                 % Names for the stores
            obj.FluxNames  = {"pn", "en", "ef", "ps", "es", "perc",...
                              "q9", "q1", "fr", "fq", "qr", "qt", "ex"};   % Names for the fluxes
            
            obj.FluxGroups.Ea = [3 5];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [12];                                      % Index or indices of fluxes to add to Streamflow
            obj.FluxGroups.Exchange  = -13;                                % Index or indices of exchange fluxes
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            x1 = theta(1);     % Maximum soil moisture storage [mm]
            x3 = theta(3);     % Maximum routing store storage [mm]
            x4 = theta(4);     % Flow delay [d]
            
            % max of stores
            obj.store_max = [x1, x3];
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh_q9 = uh_1_half(x4,delta_t);
            uh_q1 = uh_2_full(2*x4,delta_t);
            
            obj.uhs = {uh_q9, uh_q1};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        % gr4j as implemented here is subtantially different that the
        % original MARRMoT: there, S1 is solved, then q1 and q9 are routed,
        % then S2 is solved, sequentially. Here, S1 and S2 are solved at
        % the same time, the results therefore are different.
        % I have implemented it in this way so that I can keep it
        % consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that S2 is balanced at all
        % steps, which is not the case in the original MARRMoT version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            x1      = theta(1);     % Maximum soil moisture storage [mm]
            x2      = theta(2);     % Water exchange coefficient [mm/d]
            x3      = theta(3);     % Maximum routing store storage [mm]
            x4      = theta(4);     % Flow delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q9 = uhs{1};
            uh_q1 = uhs{2};
            
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
            flux_pn   = max(P-Ep,0);
            flux_en   = max(Ep-P,0);
            flux_ef   = P - flux_pn;
            flux_ps   = saturation_4(S1,x1,flux_pn);
            flux_es   = evap_11(S1,x1,flux_en);
            flux_perc = percolation_3(S1,x1);
            flux_q9   = route(.9.*(flux_pn - flux_ps + flux_perc), uh_q9);
            flux_q1   = route(.1.*(flux_pn - flux_ps + flux_perc), uh_q1);
            flux_fr   = recharge_2(3.5,S2,x3,x2);
            flux_fq   = flux_fr;
            flux_qr   = baseflow_3(S2,x3);
            flux_qt   = flux_qr + max(flux_q1 + flux_fq,0);
            % this flux is not included in original MARRMoT,
            % but it is useful to calculate the water balance
            flux_ex = flux_fr + max(flux_q1 + flux_fq,0) - flux_q1;      

            % stores ODEs
            dS1 = flux_ps - flux_es - flux_perc;
            dS2 = flux_q9 + flux_fr - flux_qr;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_pn,   flux_en, flux_ef, flux_ps, flux_es,...
                      flux_perc, flux_q9, flux_q1, flux_fr, flux_fq,...
                      flux_qr,   flux_qt, flux_ex];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_q9 = uhs{1};
            uh_q1 = uhs{2};
            
            % input fluxes to the unit hydrographs 
            fluxes = obj.fluxes(obj.t,:);  
            flux_pn   = fluxes(1);
            flux_ps   = fluxes(4);
            flux_perc = fluxes(6);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh_q9 = update_uh(uh_q9, .9.*(flux_pn - flux_ps + flux_perc));
            uh_q1 = update_uh(uh_q1, .1.*(flux_pn - flux_ps + flux_perc));
            
            obj.uhs = {uh_q9, uh_q1};
            
        end
    end
end