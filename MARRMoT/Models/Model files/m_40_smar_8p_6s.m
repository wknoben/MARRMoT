classdef m_40_smar_8p_6s < MARRMoT_model
% Class for hydrologic conceptual model: SMAR

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% O’Connell, P. E., Nash, J. E., & Farrell, J. P. (1970). River flow 
% forecasting through conceptual models part II - the Brosna catchment at 
% Ferbane. Journal of Hydrology, 10, 317–329.
%
% Tan, B. Q., & O’Connor, K. M. (1996). Application of an empirical 
% infiltration equation in the SMAR conceptual model. Journal of Hydrology,
% 185(1-4), 275–295. http://doi.org/10.1016/0022-1694(95)02993-1

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_40_smar_8p_6s()
            obj.numStores = 6;                                             % number of model stores
            obj.numFluxes = 20;                                            % number of model fluxes
            obj.numParams = 8;

            obj.JacobPattern  = [1,0,0,0,0,0;
                                 1,1,0,0,0,0;
                                 1,1,1,0,0,0;
                                 1,1,1,1,0,0;
                                 1,1,1,1,1,0;
                                 1,1,1,1,1,1];                             % Jacobian matrix of model store ODEs

            obj.parRanges = [   0,   1;     % h, Maximum fraction of direct runoff [-] 
                                0, 200;     % y, Infiltration rate [mm/d] 
                                1,2000;     % smax, Maximum soil moisture storage [mm]
                                0,   1;     % c, Evaporation reduction coefficient [-]
                                0,   1;     % g, Groundwater recharge coefficient [-]
                                0,   1;     % kg, Groundwater time parameter [d-1]
                                1,  10;     % n, Number of Nash cascade reservoirs [-]
                                1, 120];    % n*k, Routing delay [d]
        
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5" "S6"};              % Names for the stores
            obj.FluxNames  = {"pstar", "estar", "evap",   "r1", "i",...
                              "r2",    "e1",    "e2",     "e3", "e4",...
                              "e5",    "q1",    "q2",     "q3", "q4",...
                              "r3",    "rg",    "r3star", "qr", "qg"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = [3 7 8 9 10 11];                           % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [19 20];                                   % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            n    = theta(7);     % Number of Nash cascade reservoirs [-]
            nk   = theta(8);     % Routing delay [d]. n and k are optimized together 
                                 % due to interactions between them
            k    = nk/n;         % time parameter in the gamma function 
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_6_gamma(n,k,delta_t);
            
            obj.uhs = {uh};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            h    = theta(1);     % Maximum fraction of direct runoff [-] 
            y    = theta(2);     % Infiltration rate [mm/d] 
            smax = theta(3);     % Maximum soil moisture storage [mm]
            c    = theta(4);     % Evaporation reduction coefficient [-]
            g    = theta(5);     % Groundwater recharge coefficient [-]
            kg   = theta(6);     % Groundwater time parameter [d-1]
            
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
            S6 = S(6);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pstar  = effective_1(P,Ep);
            flux_estar  = effective_1(Ep,P);
            flux_evap   = min(Ep,P);
            flux_r1     = saturation_6(h,(S1+S2+S3+S4+S5),smax,flux_pstar);
            flux_i      = infiltration_4(flux_pstar-flux_r1,y);
            flux_r2     = effective_1(flux_pstar-flux_r1,flux_i);
            flux_e1     = evap_13(c,0,flux_estar,S1,delta_t);
            flux_e2     = evap_14(c,1,flux_estar,S2,S1,0.1,delta_t);
            flux_e3     = evap_14(c,2,flux_estar,S3,S2,0.1,delta_t);
            flux_e4     = evap_14(c,3,flux_estar,S4,S3,0.1,delta_t);
            flux_e5     = evap_14(c,4,flux_estar,S5,S4,0.1,delta_t);
            flux_q1     = saturation_1(flux_i, S1,smax/5);
            flux_q2     = saturation_1(flux_q1,S2,smax/5);
            flux_q3     = saturation_1(flux_q2,S3,smax/5);
            flux_q4     = saturation_1(flux_q3,S4,smax/5);
            flux_r3     = saturation_1(flux_q4,S5,smax/5);
            flux_rg     = split_1(g,flux_r3);
            flux_r3star = split_1(1-g,flux_r3);
            flux_qr     = route(flux_r1+flux_r2+flux_r3star, uh);
            flux_qg     = baseflow_1(kg,S6);

            % stores ODEs
            dS1 = flux_i  - flux_e1 - flux_q1;
            dS2 = flux_q1 - flux_e2 - flux_q2;    
            dS3 = flux_q2 - flux_e3 - flux_q3;
            dS4 = flux_q3 - flux_e4 - flux_q4;
            dS5 = flux_q4 - flux_e5 - flux_r3;
            dS6 = flux_rg - flux_qg;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5 dS6];
            fluxes = [flux_pstar, flux_estar, flux_evap,   flux_r1, flux_i,...
                      flux_r2,    flux_e1,    flux_e2,     flux_e3, flux_e4,...
                      flux_e5,    flux_q1,    flux_q2,     flux_q3, flux_q4,...
                      flux_r3,    flux_rg,    flux_r3star, flux_qr, flux_qg];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % input fluxes to the unit hydrographs
            fluxes = obj.fluxes(obj.t,:);
            flux_r1     = fluxes(4);
            flux_r2     = fluxes(6);
            flux_r3star = fluxes(18);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh = update_uh(uh, flux_r1+flux_r2+flux_r3star);
            obj.uhs = {uh};
        end
    end
end