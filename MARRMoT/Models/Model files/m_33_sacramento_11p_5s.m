classdef m_33_sacramento_11p_5s < MARRMoT_model
    % Class for sacramento model
    properties
        % in case the model has any specific properties (eg derived theta,
        % add it here)
        theta_derived   % derived parameters
    end
    methods
        
        % this function runs once as soon as the model object is created
        % and sets all the static properties of the model
        function obj = m_33_sacramento_11p_5s(delta_t, theta)          
            obj.numStores = 5;                                             % number of model stores
            obj.numFluxes = 20;                                            % number of model fluxes
            obj.numParams = 11;                                            % number of model parameters
            
            obj.JacobPattern = [1,1,0,0,0;
                                1,1,1,1,1;
                                1,1,1,1,1;
                                0,1,1,1,1;
                                0,1,1,1,1];                                % Pattern of the Jacobian matrix of model store ODEs
                            
            obj.parRanges = [0       , 1;        % pctim, Fraction impervious area [-]
                             1       , 2000;     % smax, Maximum total storage depth [mm]
                             0.005   , 0.995;    % f1, fraction of smax that is Maximum upper zone tension water storage [mm]
                             0.005   , 0.995;    % f2, fraction of smax-S1max that is Maximum upper zone free water storage [mm]
                             0       , 1;        % kuz, Interflow runoff coefficient [d-1]
                             0       , 7;        % rexp, Base percolation rate non-linearity factor [-]
                             0.005   , 0.995;    % f3, fraction of smax-S1max-S2max that is  Maximum lower zone tension water storage [mm]
                             0.005   , 0.995;    % f4, fraction of smax-S1max-S2max-S3max that is  Maximum lower zone primary free water storage [mm]
                             0       , 1;        % pfree, Fraction of percolation directed to free water stores [-]
                             0       , 1;        % klzp, Primary baseflow runoff coefficient [d-1]
                             0       , 1];       % klzs, Supplemental baseflow runoff coefficient [d-1]
            
            obj.StoreNames = ["S1", "S2", "S3", "S4", "S5"];               % Names for the stores
            obj.FluxNames  = ["qdir", "peff", "ru", "euztw", "twexu",...
                              "qsur", "qint", "euzfw", "pc", "pctw",...
                              "elztw", "twexl", "twexlp", "twexls", "pcfwp",...
                              "pcfws", "rlp", "rls", "qbfp", "qbfs"];      % Names for the fluxes
            
            obj.Flux_Ea_idx = [4 8 11];                                    % Index or indices of fluxes to add to Actual ET
            obj.Flux_Q_idx  = [1 6 7 19 20];                               % Index or indices of fluxes to add to Streamflow
            
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
            theta = obj.theta;
            smax    = theta(2);     % Maximum total storage depth [mm]
            f1      = theta(3);     % fraction of smax that is Maximum upper zone tension water storage (uztwm) [-]
            f2      = theta(4);     % fraction of smax-uztwm that is Maximum upper zone free water storage (uzfwm) [-]
            f3      = theta(7);     % fraction of smax-uztwm-uzfwm that is Maximum lower zone tension water storage (lztwm) [-]
            f4      = theta(8);     % fraction of smax-uztwm-uzfwm-lztwm that is Maximum lower zone primary free water storage (lzfwpm) [-]
            klzp    = theta(10);    % Primary baseflow runoff coefficient [d-1]
            klzs    = theta(11);    % Supplemental baseflow runoff coefficient [d-1]
            
            % Derived parameters
            % Note: we need to include some lower boundaries for store values. If
            % stores are allowed to vary freely within parameter ranges, it is possible
            % to get very small (1E-10) lower stores. The root-finding methods break
            % down in these cases, hence we need to enforce certain minimum store sizes.
            % In the most extreme case, uztwm = 0.995*1, giving 0.005/4 as mininimum 
            % size for the other stores (assuming Smax = [1,2000]).
            uztwm   = f1*smax;                                             % Maximum upper zone tension water storage [mm]
            uzfwm   = max(0.005/4,f2*(smax-uztwm));                        % Maximum upper zone free water storage [mm]
            lztwm   = max(0.005/4,f3*(smax-uztwm-uzfwm));                  % Maximum lower zone tension water storage [mm]
            lzfwpm  = max(0.005/4,f4*(smax-uztwm-uzfwm-lztwm));            % Maximum lower zone primary free water storage [mm]
            lzfwsm  = max(0.005/4,(1-f4)*(smax-uztwm-uzfwm-lztwm));        % Maximum lower zone supplemental free water storage [mm]
            pbase   = lzfwpm*klzp + lzfwsm*klzs;                           % Base percolation rate [mm/d]
            zperc   = min(100000,...
                          (lztwm+lzfwsm*(1-klzs))/(lzfwsm*klzs+lzfwpm*klzp) + ...
                          (lzfwpm*(1-klzp))/(lzfwsm*klzs+lzfwpm*klzp));    % Base percolation rate multiplication factor [-]: can return Inf, hence the min(10000,...)
                      
            obj.theta_derived = [uztwm, uzfwm, lztwm, lzfwpm, lzfwsm,...
                                 pbase, zperc];
            
            % min and max of stores
            obj.store_min = zeros(1,obj.numStores);
            obj.store_max = [uztwm,uzfwm,lztwm,lzfwpm,lzfwsm];
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            
            % parameters
            theta   = obj.theta;
            pctim   = theta(1);     % Fraction impervious area [-]
            kuz     = theta(5);     % Interflow runoff coefficient [d-1]
            rexp    = theta(6);     % Base percolation rate non-linearity factor [-]
            pfree   = theta(9);     % Fraction of percolation directed to free water stores [-]
            klzp    = theta(10);    % Primary baseflow runoff coefficient [d-1]
            klzs    = theta(11);    % Supplemental baseflow runoff coefficient [d-1]
            
            theta_d = obj.theta_derived;
            uztwm   = theta_d(1);   % Maximum upper zone tension water storage [mm]
            uzfwm   = theta_d(2);   % Maximum upper zone free water storage [mm]
            lztwm   = theta_d(3);   % Maximum lower zone tension water storage [mm]
            lzfwpm  = theta_d(4);   % Maximum lower zone primary free water storage [mm]
            lzfwsm  = theta_d(5);   % Maximum lower zone supplemental free water storage [mm]
            pbase   = theta_d(6);   % Base percolation rate [mm/d]
            zperc   = theta_d(7);   % Base percolation rate multiplication factor [-]
            
            % delta_t
            delta_t = obj.delta_t;
            
            %stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            S4 = S(4);
            S5 = S(5);
            
            % climate_input
            climate_in = obj.input_climate;
            P  = climate_in(1);
            Ep = climate_in(2);
            
            % fluxes functions
            % Original formulation using MARRMoT fluxes is very slow on sacramento,
            % individual functions have been explicitly coded underneath.
            flux_qdir    = split_1(pctim,P);
            flux_peff    = split_1(1-pctim,P);
            flux_ru      = soilmoisture_1(S1,uztwm,S2,uzfwm);
            flux_euztw   = evap_7(S1,uztwm,Ep,delta_t);
            flux_twexu   = saturation_1(flux_peff,S1,uztwm);
            flux_qsur    = saturation_1(flux_twexu,S2,uzfwm);
            flux_qint    = interflow_5(kuz,S2);
            flux_euzfw   = evap_1(S2,max(0,Ep-flux_euztw),delta_t);
            flux_pc      = percolation_4(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t);
            flux_pctw    = split_1(1-pfree,flux_pc);
            flux_elztw   = evap_7(S3,lztwm,max(0,Ep-flux_euztw-flux_euzfw),delta_t);
            flux_twexl   = saturation_1(flux_pctw,S3,lztwm);  
            flux_twexlp  = split_1(deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),flux_twexl);
            flux_twexls  = split_1(deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),flux_twexl);
            flux_pcfwp   = split_1(pfree*deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),flux_pc);
            flux_pcfws   = split_1(pfree*deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),flux_pc); 
            flux_rlp     = soilmoisture_2(S3,lztwm,S4,lzfwpm,S5,lzfwsm);
            flux_rls     = soilmoisture_2(S3,lztwm,S5,lzfwsm,S4,lzfwpm);   
            flux_qbfp    = baseflow_1(klzp,S4);
            flux_qbfs    = baseflow_1(klzs,S5);
            
            
            % stores ODEs
            dS1 = flux_peff   + flux_ru    - flux_euztw - flux_twexu;
            dS2 = flux_twexu  - flux_euzfw - flux_qsur  - flux_qint  - flux_ru - flux_pc;
            dS3 = flux_pctw   + flux_rlp   + flux_rls   - flux_elztw - flux_twexl;
            dS4 = flux_twexlp + flux_pcfwp - flux_rlp   - flux_qbfp;
            dS5 = flux_twexls + flux_pcfws - flux_rls   - flux_qbfs;
            
            % outputs
            dS = [dS1, dS2, dS3, dS4, dS5];
            fluxes = [flux_qdir, flux_peff, flux_ru, flux_euztw, flux_twexu,...
                      flux_qsur, flux_qint, flux_euzfw, flux_pc, flux_pctw,...
                      flux_elztw, flux_twexl, flux_twexlp, flux_twexls, flux_pcfwp,...
                      flux_pcfws, flux_rlp, flux_rls, flux_qbfp, flux_qbfs];
        end
        
        % STEP runs at the end of every timestep, use it to update
        % still-to-flow vectors from unit hydrographs
        function step(obj, fluxes)
            % sacramento has no unit hydrographs to update, so this is
            % empty, it still needs to exists though as it is called in the
            % MARRMoT_model.run function.
        end
    end
end