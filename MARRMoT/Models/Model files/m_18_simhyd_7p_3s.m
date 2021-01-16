classdef m_18_simhyd_7p_3s < MARRMoT_model
    % Class for simhyd model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_18_simhyd_7p_3s(delta_t, theta)          
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 10;                                            % number of model fluxes
            obj.numParams = 7;                                             % number of model parameters
            
            obj.JacobPattern  = [1,0,0;
                                 1,1,0;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [0,     5;      % INSC, Maximum interception capacity, [mm]
                             0,     600;    % COEFF, Maximum infiltration loss parameter, [mm]
                             0,     15;     % SQ, Infiltration loss exponent, [-]
                             1,     2000;   % SMSC, Maximum soil moisture capacity, [mm]
                             0,     1;      % SUB, Proportionality constant, [-]
                             0,     1;      % CRAK, Proportionality constant, [-]
                             0,     1];     % K, Slow flow time scale, [d-1]  
            
            obj.StoreNames = ["S1", "S2", "S3"];                           % Names for the stores
            obj.FluxNames  = ["Ei", "EXC", "INF", "INT", "REC",...
                              "Et", "GWF", "BAS", "SRUN", "Qt"];           % Names for the fluxes
            
            obj.Flux_Ea_idx = [1 6];                                       % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = 10;                                          % Index or indices of fluxes to add to Streamflow
            
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
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = repelem(Inf, 1, obj.numStores);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta   = obj.theta;
            insc    = theta(1);     % Maximum interception capacity, [mm]
            coeff   = theta(2);     % Maximum infiltration loss parameter, [-]
            sq      = theta(3);     % Infiltration loss exponent, [-]
            smsc    = theta(4);     % Maximum soil moisture capacity, [mm]
            sub     = theta(5);     % Proportionality constant, [-],
            crak    = theta(6);     % Proportionality constant, [-]
            k       = theta(7);     % Slow flow time scale, [d-1]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            
            % climate input
            climate_in = obj.input_climate;            
            P  = climate_in(1);
            Ep = climate_in(2);
            
            % fluxes functions
            % Original formulation using MARRMoT fluxes is slow on simhyd,
            % individual functions have been explicitly coded underneath.
%             flux_Ei   = evap_1(S1,Ep,delta_t);
%             flux_EXC  = interception_1(P,S1,insc);
%             flux_INF  = infiltration_1(coeff,sq,S2,smsc,flux_EXC);
%             flux_INT  = interflow_1(sub,S2,smsc,flux_INF);
%             flux_REC  = recharge_1(crak,S2,smsc,(flux_INF-flux_INT));
%             flux_Et   = evap_2(10,S2,smsc,Ep,delta_t);
%             flux_GWF  = saturation_1((flux_INF-flux_INT-flux_REC),S2,smsc);
%             flux_BAS  = baseflow_1(k,S3);
%             flux_SRUN = flux_EXC - flux_INF; 
%             flux_Qt   = flux_SRUN + flux_INT + flux_BAS;
            
            flux_Ei   = min(S1/delta_t, Ep);
            flux_EXC  = P.*(1-(1./(1+exp((S1-insc+.05*insc)/(.01*((.01*insc==0)+insc))))));
            flux_INF  = min(coeff.*exp((-1*sq*S2)./smsc),flux_EXC);
            flux_INT  = sub*S2/smsc*flux_INF;
            flux_REC  = crak*S2/smsc*(flux_INF-flux_INT);
            flux_Et   = min([10*S2/smsc,Ep,S2/delta_t]);
            %flux_SMF is not explicitly calculated and reported in original
            %MARRMoT, so it's not reported here either to maintain results
            %comparable
            flux_SMF  = flux_INF-flux_INT-flux_REC;
            flux_GWF  = flux_SMF.*(1-(1./(1+exp((S2-smsc+.05*smsc)/(.01*((.01*smsc==0)+smsc))))));
            flux_BAS  = k.*S3;
            flux_SRUN = flux_EXC - flux_INF; 
            flux_Qt   = flux_SRUN + flux_INT + flux_BAS;
            
            % stores ODEs
            dS1  = P        - flux_Ei  - flux_EXC;
            dS2  = flux_SMF - flux_Et  - flux_GWF;   
            dS3  = flux_REC + flux_GWF - flux_BAS;
            
            % outputs
            dS = [dS1, dS2, dS3];
            fluxes = [flux_Ei, flux_EXC, flux_INF, flux_INT, flux_REC,...
                     flux_Et, flux_GWF, flux_BAS, flux_SRUN, flux_Qt];
        end
        
        function step(obj, fluxes)
            % simhyd has no unit hydrographs to update, so this is
            % empty, it still needs to exists though as it is called in the
            % MARRMoT_model.run function.
        end
    end
end
