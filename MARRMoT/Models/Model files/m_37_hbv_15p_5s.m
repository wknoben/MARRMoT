classdef m_37_hbv_15p_5s < MARRMoT_model
    % Class for hbv model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_37_hbv_15p_5s(delta_t, theta)
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 13;                                            % number of model fluxes
            obj.numParams = 15;
            
            obj.JacobPattern  = [1,1,0,0,0;
                                 1,1,0,0,0;
                                 1,1,1,1,0;
                                 1,1,1,1,0;
                                 0,0,0,1,1];                               % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [-3 ,   5;      % TT, threshold temperature for snowfall [oC] 
                              0 ,   17;     % TTI, interval length of rain-snow spectrum [oC]
                             -3 ,   3;      % TTM, threshold temperature for snowmelt [oC]
                              0 ,   1;      % CFR, coefficient of refreezing of melted snow [-]
                              0,   20;      % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
                              0 ,   1;      % WHC, maximum water holding content of snow pack [-]
                              0 ,   4;      % CFLUX, maximum rate of capillary rise [mm/d]
                              1 ,   2000;   % FC, maximum soil moisture storage [mm]
                              0.05,    0.95;% LP, wilting point as fraction of FC [-]
                              0 ,   10;     % BETA, non-linearity coefficient of upper zone recharge [-]
                              0 ,   1;      % K0, runoff coefficient from upper zone [d-1] 
                              0 ,   4;      % ALPHA, non-linearity coefficient of runoff from upper zone [-]
                              0 ,   20;     % PERC, maximum rate of percolation to lower zone [mm/d]
                              0 ,   1;      % K1, runoff coefficient from lower zone [d-1]
                              1 ,   120];   % MAXBAS, flow routing delay [d]                 
            
            obj.StoreNames = ["S1" "S2" "S3" "S4" "S5"];                   % Names for the stores
            obj.FluxNames  = ["sf", "refr", "melt", "rf",   "in", "se", "cf",...
                              "ea", "r",    "q0",   "perc", "q1", "qt"];   % Names for the fluxes
            
            obj.Flux_Ea_idx = 8;                                           % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = 13;                                          % Index or indices of fluxes to add to Streamflow
            
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
            %parameters
            theta   = obj.theta;
            delta_t = obj.delta_t;
            maxbas  = theta(15);    % MAXBAS, flow routing delay [d]
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = inf(1,obj.numStores);
            
            % initialise the unit hydrographs and still-to-flow vectors            
            uh = uh_4_full(maxbas,delta_t);
            
            obj.uhs        = {uh};
            obj.fluxes_stf = arrayfun(@(n) zeros(1, n), cellfun(@length, obj.uhs), 'UniformOutput', false);
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation        
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            tt      = theta(1);     % TT, middle of snow-rain interval [oC]
            tti     = theta(2);     % TTI, interval length of rain-snow spectrum [oC]
            ttm     = theta(3);     % TTM, threshold temperature for snowmelt [oC]
            cfr     = theta(4);     % CFR, coefficient of refreezing of melted snow [-]
            cfmax   = theta(5);     % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
            whc     = theta(6);     % WHC, maximum water holding content of snow pack [-]
            cflux   = theta(7);   	% CFLUX, maximum rate of capillary rise [mm/d]
            fc      = theta(8);     % FC, maximum soil moisture storage [mm]
            lp      = theta(9);     % LP, wilting point as fraction of FC [-]
            beta    = theta(10);    % BETA, non-linearity coefficient of upper zone recharge [-]
            k0      = theta(11);    % K0, runoff coefficient from upper zone [d-1], 
            alpha   = theta(12);    % ALPHA, non-linearity coefficient of runoff from upper zone [-]
            perc    = theta(13);    % PERC, maximum rate of percolation to lower zone [mm/d]
            k1      = theta(14);    % K1, runoff coefficient from lower zone [d-1]
            maxbas  = theta(15);    % MAXBAS, flow routing delay [d]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % climate input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
%             flux_sf   = snowfall_2(P,T,tt,tti);
%             flux_refr = refreeze_1(cfr,cfmax,ttm,T,S2,delta_t);
%             flux_melt = melt_1(cfmax,ttm,T,S1,delta_t);
%             flux_rf   = rainfall_2(P,T,tt,tti);
%             flux_in   = infiltration_3(flux_rf+flux_melt,S2,whc*S1);
%             flux_se   = excess_1(obj.Sold(2),whc*S1,delta_t);
%             flux_cf   = capillary_1(cflux,S3,fc,S4,delta_t);
%             flux_ea   = evap_3(lp,S3,fc,Ep,delta_t);
%             flux_r    = recharge_2(beta,S3,fc,flux_in+flux_se);
%             flux_q0   = interflow_2(k0,S4,alpha,delta_t);
%             flux_perc = percolation_1(perc,S4,delta_t);
%             flux_q1   = baseflow_1(k1,S5);
%             flux_qt   = uh(1) * (flux_q0 + flux_q1) + stf(1);

            flux_sf   = min(P,max(0,P.*(tt+0.5*tti-T)/tti));
            flux_refr = min(max(0,cfr*cfmax*(ttm-T)),S2/delta_t);
            flux_melt = min(max(cfmax*(T-ttm),0),S1/delta_t);
            flux_rf   = min(P,max(0,P.*(T-(tt-0.5*tti))/tti));
            flux_in   = (flux_rf+flux_melt).*(1-...
                (1./(1+exp((S2-whc*S1+.05*whc*S1)/(.01*((.01*whc*S1==0)+(whc*S1)))))));
            flux_se   = max((obj.Sold(2)-whc*S1)/delta_t,0);
            flux_cf   = min(cflux.*(1-S3/fc),S4/delta_t);
            flux_ea   = min([S3/(lp*fc)*Ep,Ep,S3/delta_t]);
            flux_r    = (flux_in+flux_se)*((max(S3,0)/fc)^beta);
            flux_q0   = min(k0*max(S4,0)^(1+alpha),max(S4/delta_t,0));
            flux_perc = min(perc,S4/delta_t);
            flux_q1   = k1.*S5;
            flux_qt   = uh(1) * (flux_q0 + flux_q1) + stf(1);

            % stores ODEs
            dS1 = flux_sf   + flux_refr - flux_melt;
            dS2 = flux_rf   + flux_melt - flux_refr - flux_in - flux_se;
            dS3 = flux_in   + flux_se   + flux_cf   - flux_ea - flux_r;
            dS4 = flux_r    - flux_cf   - flux_q0   - flux_perc;
            dS5 = flux_perc - flux_q1;
            
            % outputs
            dS = [dS1 dS2 dS3 dS4 dS5];
            fluxes = [flux_sf,   flux_refr, flux_melt, flux_rf, flux_in,...
                      flux_se,   flux_cf,   flux_ea,   flux_r,  flux_q0,...
                      flux_perc, flux_q1, flux_qt];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % unit hydrographs and still-to-flow vectors
            uhs = obj.uhs; stf = obj.fluxes_stf;
            uh = uhs{1}; stf = stf{1};
            
            % input fluxes to the unit hydrographs 
            flux_q0 = fluxes(10);
            flux_q1 = fluxes(12);
            
            % update still-to-flow vectors using fluxes at current step and
            % unit hydrographs
            stf      = (uh .* (flux_q0 + flux_q1)) + stf;
            stf      = circshift(stf,-1);
            stf(end) = 0;
            
            obj.fluxes_stf = {stf};
        end
    end
end