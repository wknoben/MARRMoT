classdef m_21_flexb_9p_3s < MARRMoT_model
% Class for hydrologic conceptual model: Flex-B

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

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_21_flexb_9p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 9;                                             % number of model fluxes
            obj.numParams = 9; 

            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,0,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [1,2000;        % URmax, Maximum soil moisture storage [mm]
                             0, 10;         % beta, Unsaturated zone shape parameter [-]
                             0, 1;          % D, Fast/slow runoff distribution parameter [-]
                             0, 20;         % PERCmax, Maximum percolation rate [mm/d]
                             0.05, 0.95;    % Lp, Wilting point as fraction of s1max [-]
                             1, 5;          % Nlagf, Flow delay before fast runoff [d]
                             1, 15;         % Nlags, Flow delay before slow runoff [d]
                             0, 1;          % Kf, Fast runoff coefficient [d-1]
                             0, 1];         % Ks, Slow runoff coefficient [d-1]
            
            obj.StoreNames = {"S1", "S2" "S3"};                             % Names for the stores
            obj.FluxNames  = {"ru", "eur",  "ps", "rf", "rs"...
                              "rfl", "rsl", "qf", "qs"};                   % Names for the fluxes
            
            obj.FluxGroups.Ea = 2;                                         % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9];                                     % Index or indices of fluxes to add to Streamflow
            
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
        % flexb as implemented here is subtantially different that the
        % original MARRMoT: there, S1 is solved, then rsl and rfl are 
        % routed, then S2 and S3 are solved, sequentially. Here, S1, S2 and
        % S3 are all solved at the same time, the results therefore are
        % different. I have implemented it in this way so that I can keep
        % it consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that S2 and S3 are balanced at
        % all steps, which is not the case in the original MARRMoT version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            s1max   = theta(1);     % Maximum soil moisture storage [mm]
            beta    = theta(2);     % Unsaturated zone shape parameter [-]
            d       = theta(3);     % Fast/slow runoff distribution parameter [-]
            percmax = theta(4);     % Maximum percolation rate [mm/d]
            lp      = theta(5);     % Wilting point as fraction of s1max [-]
            kf      = theta(8);     % Fast runoff coefficient [d-1]
            ks      = theta(9);     % Slow runoff coefficient [d-1]
            
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
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ru  = saturation_3(S1,s1max,beta,P);
            flux_eur = evap_3(lp,S1,s1max,Ep,delta_t);
            flux_ps  = percolation_2(percmax,S1,s1max,delta_t);
            flux_rf  = split_1(1-d,P-flux_ru);
            flux_rs  = split_1(d,P-flux_ru);
            flux_rfl = route(flux_rf, uh_f);
            flux_rsl = route(flux_ps + flux_rs, uh_s);
            flux_qf  = baseflow_1(kf,S2);
            flux_qs  = baseflow_1(ks,S3);

            % stores ODEs
            dS1 = flux_ru - flux_eur - flux_ps;
            dS2 = flux_rfl - flux_qf;
            dS3 = flux_rsl - flux_qs;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ru, flux_eur,  flux_ps, flux_rf, flux_rs...
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
            flux_ps = fluxes(3);
            flux_rf = fluxes(4);
            flux_rs = fluxes(5);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh_f = update_uh(uh_f, flux_rf);
            uh_s = update_uh(uh_s, flux_ps + flux_rs);
            
            obj.uhs = {uh_f, uh_s};
        end
    end
end