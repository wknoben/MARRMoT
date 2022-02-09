classdef m_15_plateau_8p_2s < MARRMoT_model
% Class for hydrologic conceptual model: Plateau model

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Savenije, H. H. G. (2010). “Topography driven conceptual modelling 
% (FLEX-Topo).” Hydrology and Earth System Sciences, 14(12), 2681–2692. 
% https://doi.org/10.5194/hess-14-2681-2010

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_15_plateau_8p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 9;                                             % number of model fluxes
            obj.numParams = 8; 

            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0   , 200;    % Fmax, maximum infiltration rate [mm/d]
                             0   , 5;      % Dp, interception capacity [mm]
                             1   , 2000;   % SUmax, soil misture depth [mm]
                             0.05,    0.95;% Swp, wilting point as fraction of Sumax [-]
                             0   , 1;      % p, coefficient for moisture constrained evaporation [-]
                             1   , 120;    % tp, time delay for routing [d]
                             0   , 4;      % c, capillary rise [mm/d]
                             0   , 1];     % kp, base flow time parameter [d-1]
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"pe", "ei", "pie", "pi",...
                              "et", "r",  "c",   "qpgw", "qpieo"};         % Names for the fluxes
            
            obj.FluxGroups.Ea = [2 5];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9];                                     % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            tp    = theta(6);     % Time delay of surface flow [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_3_half(tp,delta_t);
            
            obj.uhs = {uh};
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            fmax  = theta(1);     % Maximum infiltration rate [mm/d]
            dp    = theta(2);     % Daily interception [mm]
            sumax = theta(3);     % Maximum soil moisture storage [mm]
            lp    = theta(4);     % Wilting point [-], defined as lp*sumax 
            p     = theta(5);     % Parameter for moisture constrained evaporation [-]
            tp    = theta(6);     % Time delay of surface flow [d]
            c     = theta(7);     % Rate of capillary rise [mm/d]
            kp    = theta(8);     % Groundwater runoff coefficient [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
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
            flux_pe    = interception_2(P,dp);
            flux_ei    = P - flux_pe;                                      % track 'intercepted' water
            flux_pi    = infiltration_4(flux_pe,fmax);
            flux_pie   = flux_pe - flux_pi;
            flux_et    = evap_4(Ep,p,S1,lp,sumax,delta_t);
            flux_c     = capillary_2(c,S2,delta_t);
            flux_r     = saturation_1((flux_pi+flux_c),S1,sumax);
            flux_qpgw  = baseflow_1(kp,S2);
            flux_qpieo = route(flux_pie, uh);

            % stores ODEs
            dS1 = flux_pi + flux_c - flux_et - flux_r;
            dS2 = flux_r - flux_c - flux_qpgw;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_pe, flux_ei, flux_pie, flux_pi,...
                      flux_et, flux_r,  flux_c,   flux_qpgw, flux_qpieo];
        end
        
        % STEP runs at the end of every timestep.
        function obj = step(obj)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs;
            uh = uhs{1};
            
            % input fluxes to the unit hydrographs 
            fluxes = obj.fluxes(obj.t,:); 
            flux_pie = fluxes(3);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            uh = update_uh(uh, flux_pie);
            obj.uhs = {uh};
            
        end
    end
end