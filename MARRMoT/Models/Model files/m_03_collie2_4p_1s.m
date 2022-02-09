classdef m_03_collie2_4p_1s < MARRMoT_model
% Class for hydrologic conceptual model: Collie River v2

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Jothityangkoon, C., M. Sivapalan, and D. Farmer (2001), Process controls
% of water balance variability in a large semi-arid catchment: downward 
% approach to hydrological model development. Journal of Hydrology, 254,
% 174198. doi: 10.1016/S0022-1694(01)00496-6.

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_03_collie2_4p_1s()
            obj.numStores = 1;                                             % number of model stores
            obj.numFluxes = 4;                                             % number of model fluxes
            obj.numParams = 4; 

            obj.JacobPattern  = [1];                                       % Jacobian matrix of model store ODEs
            
            obj.parRanges = [1   , 2000;      % Smax [mm]
                             0.05, 0.95;      % fc as fraction of Smax
                             0   , 1 ;        % a, subsurface runoff coefficient [d-1]
                             0.05, 0.95];     % M, fraction forest cover [-]
            
            obj.StoreNames = {"S1"};                                       % Names for the stores
            obj.FluxNames  = {"eb", "ev", "qse", "qss"};                   % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 2];                                     % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [3 4];                                     % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            S1max   = theta(1);     % Maximum soil moisture storage [mm]
            Sfc     = theta(2);     % Field capacity as fraction of S1max [-] 
            a       = theta(3);     % Subsurface runoff coefficient [d-1]
            M       = theta(4);     % Fraction forest cover [-]
            
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
            flux_eb   = evap_7(S1,S1max,(1-M)*Ep,delta_t);
            flux_ev   = evap_3(Sfc,S1,S1max,M*Ep,delta_t);
            flux_qse  = saturation_1(P,S1,S1max);
            flux_qss  = interflow_8(S1,a,Sfc*S1max);

            % stores ODEs
            dS1 = P - flux_eb - flux_ev - flux_qse - flux_qss;
            
            % outputs
            dS = [dS1];
            fluxes = [flux_eb, flux_ev, flux_qse, flux_qss];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end