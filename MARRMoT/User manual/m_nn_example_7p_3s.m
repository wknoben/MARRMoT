classdef m_nn_example_7p_3s < MARRMoT_model
    % Class for MARRMoT model
    
% Hydrologic conceptual model: [MARRMoT User Manual example model] 
%   
% Model reference
% MARRMoT User Manual, 2018. 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
    
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_nn_example_7p_3s(delta_t, theta)
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 9;                                             % number of model fluxes
            obj.numParams = 7;
            
            obj.JacobPattern  = [1,1,0;
                                 1,1,0;
                                 0,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [ 0,    4;     % crate, Maximum capillary rise rate [mm/d]
                              1, 2000;     % uzamx, Maximum upper zone storage [mm]
                              0,   20;     % prate, Maximum percolation rate [mm/d]
                              0,    1;     % klz, Lower zone runoff coefficient [d-1]
                              0,    1;     % alpha, Fraction of lower zone runoff to groundwater [-]
                              0,    1;     % kg, Groundwater runoff coefficient [d-1]
                              1,  120];    % d, Routing delay [d]
            
            obj.StoreNames = ["S1" "S2", "S3"];                                  % Names for the stores
            obj.FluxNames  = ["qse", "e",  "qp", "qc", "qlz",...
                              "qf",  "qg", "qs", "qt"];                    % Names for the fluxes
            
            obj.Flux_Ea_idx = 2;                                           % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = 9;                                           % Index or indices of fluxes to add to Streamflow
            
            % setting delta_t and theta triggers the function obj.init()
            if nargin > 0 && ~isempty(delta_t)
                obj.delta_t = delta_t;
            end
            if nargin > 1 && ~isempty(theta)
                obj.theta = theta;
            end
        end
        
        % INIT is run automatically as soon as both theta and delta_t are
        % set (it is therefore ran only once at the beginning of the run. 
        % Use it to initialise all the model parameters (in case there are
        % derived parameters) and unit hydrographs and set minima and
        % maxima for stores based on parameters.
        function obj = init(obj)
            theta   = obj.theta;
            delta_t = obj.delta_t;
            
            d = theta(7);     % Routing delay [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(d,delta_t);
            
            obj.uhs        = {uh};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            crate = theta(1);     % Maximum capillary rise rate [mm/d]
            uzmax = theta(2);     % Maximum upper zone storage [mm]
            prate = theta(3);     % Maximum percolation rate [mm/d]
            klz   = theta(4);     % Lower zone runoff coefficient [d-1]
            alpha = theta(5);     % Fraction of lower zone runoff to groundwater [-]
            kg    = theta(6);     % Groundwater runoff coefficient [d-1]
            d     = theta(7);     % Routing delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            
            % climate input
            c = obj.input_climate;
            P  = c(1);
            Ep = c(2);
            T  = c(3);
            
            % fluxes functions
            flux_qse = saturation_1(P,S1,uzmax);
            flux_e   = evap_7(S1,uzmax,Ep,delta_t);
            flux_qp  = percolation_1(prate,S1,delta_t);
            flux_qc  = capillary_1(crate,S1,uzmax,S2,delta_t);
            flux_qlz = baseflow_1(klz,S2);
            flux_qf  = split_1(1-alpha,flux_qlz);
            flux_qg  = split_1(alpha,flux_qlz);
            flux_qs  = baseflow_1(kg,S3);
            flux_qt  = uh(1) .* (flux_qse + flux_qf + flux_qs) + stf(1);

            % stores ODEs
            dS1 = P         + flux_qc - flux_e  - flux_qse - flux_qp;
            dS2 = flux_qp - flux_qc  - flux_qlz;
            dS3 = flux_qg - flux_qs;
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_qse, flux_e,  flux_qp, flux_qc, flux_qlz,...
                      flux_qf,  flux_qg, flux_qs, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
                        
            % input fluxes to the unit hydrographs 
            flux_qse = fluxes(1);
            flux_qf  = fluxes(6);
            flux_qs  = fluxes(8);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf      = (uh .* (flux_qse + flux_qf + flux_qs)) + stf;
            stf      = circshift(stf,-1);
            stf(end) = 0;
            
            obj.fluxes_stf = {stf};
        end
    end
end