classdef m_39_mcrm_16p_5s < MARRMoT_model
% Class for hydrologic conceptual model: Midland Catchment Runoff Model

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model references
% Moore, R. J., & Bell, V. A. (2001). Comparison of rainfall-runoff models 
% for flood forecasting. Part 1: Literature review of models. Bristol: 
% Environment Agency.

    properties
        % model-specific attributes
        
        aux_theta      % Auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_39_mcrm_16p_5s()
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 12;                                            % number of model fluxes
            obj.numParams = 16;
            
            obj.JacobPattern  = [1,0,0,0,0;
                                 1,1,0,0,0;
                                 0,1,1,0,0;
                                 1,1,1,1,0;
                                 1,1,1,1,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0   , 5;       % smax, Maximum interception storage [mm]
                             0.01, 0.99;    % cmax, Maximum fraction of area contributing to rapid runoff [-]
                             0.01, 0.99;    % ct, Fraction of cmax that is the minimum contributing area [-]
                             0   , 2;       % c1, Shape parameter for rapid flow distribution [mm-1]
                             0   , 1;       % ce, Shape parameter for evaporation [mm-1]
                             1   , 2000;    % dsurp, Threshold for direct runoff [mm]
                             0   , 1;       % kd, Direct runoff time parameter [d-1]
                             1   , 5;       % gamd, Direct runoff flow non-linearity [-]
                             0   , 20;      % qpmax, Maximum percolation rate [mm/d]
                             0   , 1;       % kg, Groundwater time parameter [d-1]
                             1   , 120;     % tau, Routing delay [d]
                             1   , 300;     % sbf, Maximum routing store depth [mm]
                             0   , 1;       % kcr, Channel flow time parameter [d-1]
                             1   , 5;       % gamcr, Channel flow non-linearity [-]
                             0   , 1;       % kor, Out-of-bank flow time parameter [d-1]
                             1   , 5];      % gamor, Out-of-bank flow non-linearity [-]
     
            obj.StoreNames = {"S1", "S2" "S3" "S4" "S5"};                   % Names for the stores
            obj.FluxNames  = {"ec", "qt", "qr",  "er",  "qn", "qd",...
                              "qp", "qb", "uib", "uob", "qic", "qoc"};     % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 4];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [11 12];                                   % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            cmax   = theta(2);      % Maximum fraction of area contributing to rapid runoff [-]
            ct     = theta(3);      % Fraction of cmax that is the minimum contributing area c0 [-]
            tau    = theta(11);     % Routing delay [d]
            
            % auxiliary parameters
            c0     = ct*cmax;       % Minimum fraction of area contributing to rapid runoff [-]
            obj.aux_theta(1) = c0;
            
            % min and max of stores
            obj.store_min = [0,-1E6,0,0,0];
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_7_uniform(tau,delta_t);
            
            obj.uhs = {uh};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        % mcrm as implemented here is subtantially different that the
        % original MARRMoT: there, S1, S2 and S3 are solved, then uib is
        % routed, then S4 and S5 are solved, sequentially. Here, all stores
        % are solved at the same time, the results therefore are different.
        % I have implemented it in this way so that I can keep it
        % consistent with other models and use a single call to
        % MARRMoT_model.solve_stores to solve the stores' ODEs, this
        % implementation actually guarantees that all stores are balanced
        % at all steps, which is not the case in the original MARRMoT
        % version.
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            smax   = theta(1);      % Maximum interception storage [mm]
            cmax   = theta(2);      % Maximum fraction of area contributing to rapid runoff [-]
            c1     = theta(4);      % Shape parameter for rapid flow distribution [mm-1]
            ce     = theta(5);      % Shape parameter for evaporation [mm-1]
            dsurp  = theta(6);      % Threshold for direct runoff [mm]
            kd     = theta(7);      % Direct runoff time parameter [d-1]
            gamd   = theta(8);      % Direct runoff flow non-linearity [-]
            qpmax  = theta(9);      % Maximum percolation rate [mm/d]
            kg     = theta(10);     % Groundwater time parameter [d-1]
            sbf    = theta(12);     % Maximum routing store depth [mm]
            kcr    = theta(13);     % Channel flow time parameter [d-1]
            gamcr  = theta(14);     % Channel flow non-linearity [-]
            kor    = theta(15);     % Out-of-bank flow time parameter [d-1]
            gamor  = theta(16);     % Out-of-bank flow non-linearity [-]
            
            % auxiliary parameters
            c0 = obj.aux_theta(1);
            
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
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_ec  = evap_1(S1,Ep,delta_t);
            flux_qt  = interception_1(P,S1,smax);
            flux_qr  = saturation_10(cmax,c0,c1,S2,flux_qt);
            flux_er  = evap_17(ce,S2,Ep-flux_ec);
            flux_qn  = effective_1(flux_qt,flux_qr);
            flux_qd  = interflow_9(S2,kd,dsurp,gamd,delta_t);
            flux_qp  = percolation_6(qpmax,dsurp,S2,delta_t);
            flux_qb  = baseflow_7(kg,1.5,S3,delta_t);
            flux_uib = route(flux_qr + flux_qd + flux_qb, uh);
            flux_uob = saturation_1(flux_uib,S4,sbf);
            flux_qic = routing_1(kcr,gamcr,3/4,S4,delta_t);
            flux_qoc = routing_1(kor,gamor,3/4,S5,delta_t);

            % stores ODEs
            dS1 = P - flux_ec - flux_qt;
            dS2 = flux_qn - flux_er - flux_qd - flux_qp;
            dS3 = flux_qp - flux_qb;
            dS4 = flux_uib - flux_uob - flux_qic;
            dS5 = flux_uob - flux_qoc;
             
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_ec,  flux_qt,  flux_qr,  flux_er,...
                      flux_qn,  flux_qd,  flux_qp,  flux_qb,...
                      flux_uib, flux_uob, flux_qic, flux_qoc];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function obj = step(obj)
           % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % input fluxes to the unit hydrographs
            fluxes = obj.fluxes(obj.t,:);
            flux_qr = fluxes(3);
            flux_qd = fluxes(6);
            flux_qb = fluxes(8);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh = update_uh(uh, flux_qr + flux_qd + flux_qb);
            obj.uhs = {uh};
        end
    end
end