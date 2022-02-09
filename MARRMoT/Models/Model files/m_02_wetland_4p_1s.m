classdef m_02_wetland_4p_1s < MARRMoT_model
% Class for hydrologic conceptual model: Wetland model

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
        function obj = m_02_wetland_4p_1s()
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 5;                                             % number of model fluxes
            obj.numParams = 4;

            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [0   , 5;      % Dw, interception capacity [mm]
                             0   , 10;     % Betaw, soil misture distribution parameter [-]
                             1   , 2000;   % Swmax, soil misture depth [mm]
                             0   , 1];     % kw, base flow time parameter [d-1]
            
            obj.StoreNames = {"S1"};                                       % Names for the stores
            obj.FluxNames  = {"pe",  "ei",  "ew",  "qwsof", "qwgw"};       % Names for the fluxes
            
            obj.FluxGroups.Ea = [2 3];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [4 5];                                     % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            dw    = theta(1);     % Daily interception [mm]
            betaw = theta(2);     % Soil moisture storage distribution parameter [-]
            swmax = theta(3);     % Maximum soil moisture storage [mm], 
            kw    = theta(4);     % Runoff coefficient [d-1]
            
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            flux_pe    = interception_2(P,dw);
            flux_ei    = P - flux_pe;
            flux_ew    = evap_1(S1,Ep,delta_t);
            flux_qwsof = saturation_2(S1,swmax,betaw,flux_pe);
            flux_qwgw  = baseflow_1(kw,S1);

            % stores ODEs
            dS1 = flux_pe - flux_ew - flux_qwsof - flux_qwgw; 
            
            % outputs
            dS = [dS1];
            fluxes = [flux_pe,  flux_ei,  flux_ew,  flux_qwsof, flux_qwgw];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end