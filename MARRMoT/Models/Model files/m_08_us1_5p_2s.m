classdef m_08_us1_5p_2s < MARRMoT_model
% Class for hydrologic conceptual model: United States model v1

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Bai, Y., Wagener, T., & Reed, P. (2009). A top-down framework for 
% watershed model evaluation and selection under uncertainty. Environmental
% Modelling & Software, 24(8), 901–916. 
% http://doi.org/10.1016/j.envsoft.2008.12.012

    properties
        % model-specific attributes
    end
    methods
        
        % creator method
        function obj = m_08_us1_5p_2s()
            obj.numStores = 2;                                             % number of model stores
            obj.numFluxes = 9;                                             % number of model fluxes
            obj.numParams = 5; 

            obj.JacobPattern  = [1,1;
                                 1,1];                                     % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0   , 1;       % Alpha_ei, Fraction of intercepted rainfall [-]
                             0.05, 0.95;    % M, Fraction forest [-]
                             1   , 2000;    % Smax, Maximum soil moisture [mm]
                             0.05, 0.95;    % fc, Field capacity as fraction of Smax [-]
                             0   , 1];      % Alpha_ss, Subsurface routing delay [d-1]
            
            obj.StoreNames = {"S1", "S2"};                                  % Names for the stores
            obj.FluxNames  = {"eusei",  "eusveg", "eusbs", "esatveg",...
                              "esatbs", "rg",     "se",    "qse", "qss"};  % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 2 3 4 5];                               % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [8 9];                                     % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            alpha_ei    = theta(1);     % Fraction of intercepted rainfall [-]
            m           = theta(2);     % Fraction forest [-]
            smax        = theta(3);     % Maximum soil moisture [mm]
            fc          = theta(4);     % Field capacity as fraction of smax [-]
            alpha_ss    = theta(5);     % Subsurface routing delay [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
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
            flux_eusei   = interception_3(alpha_ei,P);
            flux_eusveg  = evap_8(S1,S2,m,fc*(smax-S2),Ep,delta_t);
            flux_eusbs   = evap_9(S1,S2,m,smax,Ep,delta_t);
            flux_esatveg = evap_10(m,S2,S1+S2,Ep,delta_t);
            flux_esatbs  = evap_5(m,S2,S1+S2,Ep,delta_t);
            flux_rg      = saturation_1(P,S1,fc*(smax-S2));
            flux_se      = excess_1(S1,fc*(smax-S2),delta_t);
            flux_qse     = saturation_1(flux_rg+flux_se,S2,smax);
            flux_qss     = baseflow_1(alpha_ss,S2);

            % stores ODEs
            dS1 =   P     - flux_eusei - flux_eusveg  - flux_eusbs  - flux_rg  - flux_se;
            dS2 = flux_rg + flux_se    - flux_esatveg - flux_esatbs - flux_qse - flux_qss;
            
            % outputs
            dS = [dS1 dS2];
            fluxes = [flux_eusei,   flux_eusveg, flux_eusbs, ...
                      flux_esatveg, flux_esatbs, flux_rg,...
                      flux_se,      flux_qse,    flux_qss];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end