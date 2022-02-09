classdef m_34_flexis_12p_5s < MARRMoT_model
% Class for hydrologic conceptual model: Flex-IS

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Fenicia, F., McDonnell, J. J., & Savenije, H. H. G. (2008). Learning from
% model improvement: On the contribution of complementary data to process
% understanding. Water Resources Research, 44(6), 1–13. 
% http://doi.org/10.1029/2007WR006386
%
% Nijzink, R., Hutton, C., Pechlivanidis, I., Capell, R., Arheimer, B., 
% Freer, J., … Hrachowitz, M. (2016). The evolution of root zone moisture 
% capacities after land use change: a step towards predictions under 
% change? Hydrology and Earth System Sciences Discussions, 20(August), 
% 4775–4799. http://doi.org/10.5194/hess-2016-427

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_34_flexis_12p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 14;                                             % number of model fluxes
            obj.numParams = 12;
            
            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 1,1,1,0,0;
                                 1,1,1,1,0;
                                 1,1,1,0,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1,2000;        % URmax, Maximum soil moisture storage [mm]
                             0, 10;         % beta, Unsaturated zone shape parameter [-]
                             0, 1;          % D, Fast/slow runoff distribution parameter [-]
                             0, 20;         % PERCmax, Maximum percolation rate [mm/d]
                             0.05, 0.95;    % Lp, Wilting point as fraction of s1max [-]
                             1, 5;          % Nlagf, Flow delay before fast runoff [d]
                             1, 15;         % Nlags, Flow delay before slow runoff [d]
                             0, 1;          % Kf, Fast runoff coefficient [d-1]
                             0, 1;          % Ks, Slow runoff coefficient [d-1]
                             0, 5;          % Imax, Maximum interception storage [mm]
                            -3, 5;          % TT, Threshold temperature for snowfall/snowmelt [oC]
                             0, 20];        % ddf, Degree-day factor for snowmelt [mm/d/oC]
            
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"ps", "pi", "m", "peff", "ei",...
                              "ru", "eur", "rp", "rf", "rs",...
                              "rf1", "rs1", "qf", "qs"};                   % Names for the fluxes
            
            obj.FluxGroups.Ea = [5 7];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [13 14];                                   % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            nlagf   = theta(6);     % Flow delay before fast runoff [d]
            nlags   = theta(7);     % Flow delay before slow runoff [d]
            
            % initialise the unit hydrographs and still-to-flow vectors
            uh_f = uh_3_half(nlagf,delta_t);
            uh_s = uh_3_half(nlags,delta_t);
            
            obj.uhs = {uh_f, uh_s};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        % flexis as implemented here is subtantially different that the
        % original MARRMoT: there, S1, S2 and S3 are solved, then rfl and 
        % rsl are routed, then S4 and S5 are solved, sequentially. Here,
        % all stores are solved at the same time, the results therefore are
        % different. I have implemented it in this way so that I can keep
        % it consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that all stores are balanced
        % at all steps, which is not the case in the original MARRMoT
        % version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            smax    = theta(1);     % Maximum soil moisture storage [mm]
            beta    = theta(2);     % Unsaturated zone shape parameter [-]
            d       = theta(3);     % Fast/slow runoff distribution parameter [-]
            percmax = theta(4);     % Maximum percolation rate [mm/d]
            lp      = theta(5);     % Wilting point as fraction of s1max [-]
            nlagf   = theta(6);     % Flow delay before fast runoff [d]
            nlags   = theta(7);     % Flow delay before slow runoff [d]
            kf      = theta(8);     % Fast runoff coefficient [d-1]
            ks      = theta(9);     % Slow runoff coefficient [d-1]
            imax    = theta(10);    % Maximum interception storage [mm]
            tt      = theta(11);    % Threshold temperature for snowfall/snowmelt [oC]
            ddf     = theta(12);    % Degree-day factor for snowmelt [mm/d/oC]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_f = uhs{1};
            uh_s = uhs{2};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ps   = snowfall_1(P,T,tt);
            flux_pi   = rainfall_1(P,T,tt);
            flux_m    = melt_1(ddf,tt,T,S1,delta_t);
            flux_peff = interception_1( flux_m+flux_pi ,S2,imax);
            flux_ei   = evap_1(S2,Ep,delta_t);
            flux_ru   = saturation_3(S3,smax,beta,flux_peff);
            flux_eur  = evap_3(lp,S3,smax,Ep,delta_t);
            flux_rp   = percolation_2(percmax,S3,smax,delta_t);
            flux_rf   = split_1(1-d,flux_peff-flux_ru);
            flux_rs   = split_1(d,flux_peff-flux_ru);
            flux_rfl = route(flux_rf, uh_f);
            flux_rsl = route(flux_rs + flux_rp, uh_s);
            flux_qf   = baseflow_1(kf,S4);
            flux_qs   = baseflow_1(ks,S5);
            
            % stores ODEs
            dS1 = flux_ps  - flux_m;
            dS2 = flux_m   + flux_pi  - flux_peff - flux_ei;
            dS3 = flux_ru  - flux_eur - flux_rp;
            dS4 = flux_rfl - flux_qf;
            dS5 = flux_rsl - flux_qs; 

            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_ps,  flux_pi,  flux_m,  flux_peff, flux_ei,...
                      flux_ru,  flux_eur, flux_rp, flux_rf, flux_rs,... 
                      flux_rfl, flux_rsl, flux_qf, flux_qs];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh_f = uhs{1};
            uh_s = uhs{2};
            
            % input fluxes to the unit hydrographs
            fluxes = obj.fluxes(obj.t,:);
            flux_rp   = fluxes(8);
            flux_rf   = fluxes(9);
            flux_rs   = fluxes(10);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh_f = update_uh(uh_f, flux_rf);
            uh_s = update_uh(uh_s, flux_rs + flux_rp);
            
            obj.uhs = {uh_f, uh_s};
        end
    end
end