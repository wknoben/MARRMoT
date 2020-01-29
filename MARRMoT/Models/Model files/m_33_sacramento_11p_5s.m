function [ fluxOutput, fluxInternal, storeInternal, waterBalance  ] = ...
        m_33_sacramento_11p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: SACRAMENTO-SMA
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% National Weather Service (2005), II.3-SAC-SMA: Conceptualization of the 
% Sacramento Soil Moisture Accounting model.In National Weather Service 
% River Forecast System (NWSRFS) User Manual, 113
%
% Koren, V. I., Smith, M., Wang, D., & Zhang, Z. (2000). Use of soil 
% property data in the derivation of conceptual rainfall-runoff model 
% parameters. Proceedings of the 15th Conference on Hydrology, AMS, Long 
% Beach, CA, (1), 103–106.

% Steps
% --- Practical ---
% 0. Handle inputs
%
% --- Model setup ---
% 1. Set out ODE
% 2. Set out constitutive functions
% 3. Determine smoothing
% 4. Determine numerical scheme

% --- Model use ---
% 5. Solve 

% --- Practical ---
% 6. Handle outputs

%% Setup
%%INPUTS
% Time step size 
delta_t = fluxInput.delta_t;

% Data
P     = fluxInput.precip./delta_t;          % [mm/delta_t] -> [mm/d]       
Ep    = fluxInput.pet./delta_t;             % [mm/delta_t] -> [mm/d]
T     = fluxInput.temp;
t_end = length(P);

% Calibration parameters
% [name in documentation] = theta(order in which specified in parameter file)
pctim   = theta(1);     % Fraction impervious area [-]
smax    = theta(2);     % Maximum total storage depth [mm]
f1      = theta(3);     % fraction of smax that is Maximum upper zone tension water storage (uztwm) [-]
f2      = theta(4);     % fraction of smax-uztwm that is Maximum upper zone free water storage (uzfwm) [-]
kuz     = theta(5);     % Interflow runoff coefficient [d-1]
rexp    = theta(6);     % Base percolation rate non-linearity factor [-]
f3      = theta(7);     % fraction of smax-uztwm-uzfwm that is Maximum lower zone tension water storage (lztwm) [-]
f4      = theta(8);     % fraction of smax-uztwm-uzfwm-lztwm that is Maximum lower zone primary free water storage (lzfwpm) [-]
pfree   = theta(9);     % Fraction of percolation directed to free water stores [-]
klzp    = theta(10);    % Primary baseflow runoff coefficient [d-1]
klzs    = theta(11);    % Supplemental baseflow runoff coefficient [d-1]

% Derived parameters
% Note: we need to include some lower boundaries for store values. If
% stores are allowed to vary freely within parameter ranges, it is possible
% to get very small (1E-10) lower stores. The root-finding methods break
% down in these cases, hence we need to enforce certain minimum store sizes.
% In the most extreme case, uztwm = 0.995*1, giving 0.005/4 as mininimum 
% size for the other stores (assuming Smax = [1,2000]).
uztwm   = f1*smax;                                                          % Maximum upper zone tension water storage [mm]
uzfwm   = max(0.005/4,f2*(smax-uztwm));                                     % Maximum upper zone free water storage [mm]
lztwm   = max(0.005/4,f3*(smax-uztwm-uzfwm));                               % Maximum lower zone tension water storage [mm]
lzfwpm  = max(0.005/4,f4*(smax-uztwm-uzfwm-lztwm));                         % Maximum lower zone primary free water storage [mm]
lzfwsm  = max(0.005/4,(1-f4)*(smax-uztwm-uzfwm-lztwm));                     % Maximum lower zone supplemental free water storage [mm]
pbase   = lzfwpm*klzp + lzfwsm*klzs;                                        % Base percolation rate [mm/d]
zperc   = min(100000,...
              (lztwm+lzfwsm*(1-klzs))/(lzfwsm*klzs+lzfwpm*klzp) + ...
              (lzfwpm*(1-klzp))/(lzfwsm*klzs+lzfwpm*klzp));                 % Base percolation rate multiplication factor [-]: can return Inf, hence the min(10000,...)

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial UZ tension water
S20     = storeInitial(2);       % Initial UZ free water
S30     = storeInitial(3);       % Initial LZ tension water
S40     = storeInitial(4);       % Initial LZ primary free water
S50     = storeInitial(5);       % Initial LZ supplemental free water

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0];                        % lower bounds of stores
store_upp = [uztwm,uzfwm,lztwm,lzfwpm,lzfwsm];  % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1    = zeros(1,t_end);
store_S2    = zeros(1,t_end);
store_S3    = zeros(1,t_end);
store_S4    = zeros(1,t_end);
store_S5    = zeros(1,t_end);

flux_qdir   = zeros(1,t_end);
flux_peff   = zeros(1,t_end);
flux_ru     = zeros(1,t_end);
flux_euztw  = zeros(1,t_end);
flux_twexu  = zeros(1,t_end);
flux_qsur   = zeros(1,t_end);
flux_qint   = zeros(1,t_end);
flux_euzfw  = zeros(1,t_end);
flux_pc     = zeros(1,t_end);
flux_pctw   = zeros(1,t_end);
flux_elztw  = zeros(1,t_end);
flux_twexl  = zeros(1,t_end);
flux_twexlp = zeros(1,t_end);
flux_twexls = zeros(1,t_end);
flux_pcfwp  = zeros(1,t_end);
flux_pcfws  = zeros(1,t_end);
flux_rlp    = zeros(1,t_end);
flux_rls    = zeros(1,t_end);
flux_qbfp   = zeros(1,t_end);
flux_qbfs   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. UZ tension water
% S2. UZ free water
% S3. LZ tension water
% S4. LZ primary free water
% S5. LZ supplemental free water

% QDIR(pctim,P(t)): fraction of precip on impervious area
QDIR = split_1;

% PEFF(1-pctim,P(t)): fraction of precip that infiltrates
PEFF = split_1;

% RU(S1,uztwm,S2,uzfwm): water rebalance in upper zone
RU = soilmoisture_1;

% EUZTW(S1,uztwm,Ep(t),delta_t): evap from tension water
EUZTW = evap_7;

% TWEXU(PEFF(1-pctim,P(t)),S1,uztwm): outflow from tension water
TWEXU = saturation_1;

% EUZFW(S2,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t))),delta_t): evap from upper 
% zone free water, if Ep demand still exists
EUZFW = evap_1;

% QSUR(TWEXU(PEFF(1-pctim,P(t)),S1,uztwm),S2,uzfwm): overflow from upper 
% free water
QSUR = saturation_1;

% QINT(kuz,S2): interflow from free water
QINT = interflow_5;

% PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t): 
% percolation from upper free water based on demand
PC = percolation_4;

% PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)):
% percolation to tension water
PCTW = split_1;

% RLP(S3,lztwm,S4,lzfwpm,S5,lzfwsm): water rebalance in the lower zone
RLP = soilmoisture_2;

% RLS(S3,lztwm,S5,lzfwsm,S4,lzfwpm): water rebalance in the lower zone
RLS = soilmoisture_2;

% ELZTW(S3,lztwm,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t))-EUZFW(S2,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t)))),delta_t):
% evap from lower tension store if Ep demand still exists
ELZTW = evap_7;

% TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm)),S3,lztwm):
% lower zone tension water excess
TWEXL = saturation_1;

% PCFWP(pfree*deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)):
% fraction of percolation to primary free water store
PCFWP = split_1;

% TWEXLP(deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)),S3,lztwm)):
% fraction of tension water excess to primary free water
TWEXLP = split_1;

% QBFP(klzp,S4): primary baseflow
QBFP = baseflow_1;

% PCFWS(pfree*deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm)):
% fraction of percolation to supplemental free water store
PCFWS = split_1;

% TWEXLS(deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm)),S3,lztwm)):
% fraction of tension water excess to supplemental free water
TWEXLS = split_1;

% QBFS(klzs,S5): supplemental baseflow
QBFS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1,0,0,0;
                                               1,1,1,1,1;
                                               0,0,1,1,1;
                                               0,0,1,1,1;
                                               0,0,1,1,1]);                 % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,1,0,0,0;
                                                    1,1,1,1,1;
                                                    0,0,1,1,1;
                                                    0,0,1,1,1;
                                                    0,0,1,1,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    if t == 1; S4old = S40; else; S4old = store_S4(t-1); end                % store 4 at t-1
    if t == 1; S5old = S50; else; S5old = store_S5(t-1); end                % store 5 at t-1
    
    % Note: fluxes going to primary and supplemental lower zone free water
    % stores are supposed to be divided between the two according to their
    % current relative deficit (takes into account both current storage and
    % relative maximum size of each store). This computation is outsourced
    % to a dedicated function 'deficitBasedDistribution' that is able to
    % effectively deal with the no-deficit case, which results in a
    % divide-by-zero error.
                  
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5) ...
                (PEFF(1-pctim,P(t)) + ...
                 RU(S1,uztwm,S2,uzfwm) - ...
                 EUZTW(S1,uztwm,Ep(t),delta_t) - ...
                 TWEXU(PEFF(1-pctim,P(t)),S1,uztwm));                       % store 1 function with current flux values
   
    tmpf_S2 = @(S1,S2,S3,S4,S5) ...
                (TWEXU(PEFF(1-pctim,P(t)),S1,uztwm) - ...
                 EUZFW(S2,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t),delta_t)),delta_t) - ...
                 QSUR(TWEXU(PEFF(1-pctim,P(t)),S1,uztwm),S2,uzfwm) - ...
                 QINT(kuz,S2) - ...
                 PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t) - ...
                 RU(S1,uztwm,S2,uzfwm));                                    % Store 2 function
             
    tmpf_S3 = @(S1,S2,S3,S4,S5) ...
                (PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)) + ...
                 RLP(S3,lztwm,S4,lzfwpm,S5,lzfwsm) + ...
                 RLS(S3,lztwm,S5,lzfwsm,S4,lzfwpm) - ...
                 ELZTW(S3,lztwm,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t),delta_t)-EUZFW(S2,max(0,Ep(t)-EUZTW(S1,uztwm,Ep(t),delta_t)),delta_t)),delta_t) - ...
                 TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)),S3,lztwm)); % Store 3 function
             
    tmpf_S4 = @(S1,S2,S3,S4,S5) ...
                (TWEXLP(deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)),S3,lztwm)) + ...
                 PCFWP(pfree*deficitBasedDistribution(S4,lzfwpm,S5,lzfwsm),PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)) - ...
                 RLP(S3,lztwm,S4,lzfwpm,S5,lzfwsm) - ...
                 QBFP(klzp,S4));                                            % Store 4 function
    
    tmpf_S5 = @(S1,S2,S3,S4,S5) ...
                (TWEXLS(deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),TWEXL(PCTW(1-pfree,PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)),S3,lztwm)) + ...
                 PCFWS(pfree*deficitBasedDistribution(S5,lzfwsm,S4,lzfwpm),PC(pbase,zperc,rexp,max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5),lztwm+lzfwpm+lzfwsm,S2,uzfwm,delta_t)) - ...
                 RLS(S3,lztwm,S5,lzfwsm,S4,lzfwpm) - ...
                 QBFS(klzs,S5));                                            % Store 5 function
             
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5);             % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                            eq_sys(5)),...                                  % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old],...                 % storage values on previous time step
                        fsolve_options);                                    % solver options

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
    
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,tmp_resnorm,~] = rerunSolver('lsqnonlin', ...             % [tmp_sNew,tmp_resnorm,flag]
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),eq_sys(5)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old,S5old], ...% storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_qdir(t)    = QDIR(pctim,P(t));
    flux_peff(t)    = PEFF(1-pctim,P(t));
    flux_ru(t)      = RU(tmp_sFlux(1),uztwm,tmp_sFlux(2),uzfwm);
    flux_euztw(t)   = EUZTW(tmp_sFlux(1),uztwm,Ep(t),delta_t);
    flux_twexu(t)   = TWEXU(flux_peff(t),tmp_sFlux(1),uztwm);
    flux_qsur(t)    = QSUR(flux_twexu(t),tmp_sFlux(2),uzfwm);
    flux_qint(t)    = QINT(kuz,tmp_sFlux(2));
    flux_euzfw(t)   = EUZFW(tmp_sFlux(2),max(0,Ep(t)-flux_euztw(t)),delta_t);
    flux_pc(t)      = PC(pbase,zperc,rexp,max(0,lztwm-tmp_sFlux(3))+max(0,lzfwpm-tmp_sFlux(4))+max(0,lzfwsm-tmp_sFlux(5)),lztwm+lzfwpm+lzfwsm,tmp_sFlux(2),uzfwm,delta_t);
    flux_pctw(t)    = PCTW(1-pfree,flux_pc(t));
    flux_elztw(t)   = ELZTW(tmp_sFlux(3),lztwm,max(0,Ep(t)-flux_euztw(t)-flux_euzfw(t)),delta_t);
    flux_twexl(t)   = TWEXL(flux_pctw(t),tmp_sFlux(3),lztwm);  
    flux_twexlp(t)  = TWEXLP(deficitBasedDistribution(tmp_sFlux(4),lzfwpm,tmp_sFlux(5),lzfwsm),flux_twexl(t));
    flux_twexls(t)  = TWEXLS(deficitBasedDistribution(tmp_sFlux(5),lzfwsm,tmp_sFlux(4),lzfwpm),flux_twexl(t));
    flux_pcfwp(t)   = PCFWP(pfree*deficitBasedDistribution(tmp_sFlux(4),lzfwpm,tmp_sFlux(5),lzfwsm),flux_pc(t));
    flux_pcfws(t)   = PCFWS(pfree*deficitBasedDistribution(tmp_sFlux(5),lzfwsm,tmp_sFlux(4),lzfwpm),flux_pc(t)); 
    flux_rlp(t)     = RLP(tmp_sFlux(3),lztwm,tmp_sFlux(4),lzfwpm,tmp_sFlux(5),lzfwsm);
    flux_rls(t)     = RLS(tmp_sFlux(3),lztwm,tmp_sFlux(5),lzfwsm,tmp_sFlux(4),lzfwpm);   
    flux_qbfp(t)    = QBFP(klzp,tmp_sFlux(4));
    flux_qbfs(t)    = QBFS(klzs,tmp_sFlux(5));  
    
    % Update the stores
    store_S1(t) = S1old + (flux_peff(t)   + flux_ru(t)    - flux_euztw(t) - flux_twexu(t)) * delta_t;
    store_S2(t) = S2old + (flux_twexu(t)  - flux_euzfw(t) - flux_qsur(t)  - flux_qint(t)  - flux_ru(t) - flux_pc(t)) * delta_t;
    store_S3(t) = S3old + (flux_pctw(t)   + flux_rlp(t)   + flux_rls(t)   - flux_elztw(t) - flux_twexl(t)) * delta_t;
    store_S4(t) = S4old + (flux_twexlp(t) + flux_pcfwp(t) - flux_rlp(t)   - flux_qbfp(t)) * delta_t;
    store_S5(t) = S5old + (flux_twexls(t) + flux_pcfws(t) - flux_rls(t)   - flux_qbfs(t)) * delta_t;
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_euztw + flux_euzfw + flux_elztw) * delta_t;
    fluxOutput.Q      = (flux_qdir + flux_qsur + flux_qint + flux_qbfp + flux_qbfs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.qdir  = flux_qdir * delta_t;
    fluxInternal.qsur  = flux_qsur * delta_t;
    fluxInternal.qint  = flux_qint * delta_t;
    fluxInternal.qbfp  = flux_qbfp * delta_t;
    fluxInternal.qbfs  = flux_qbfs * delta_t;
    fluxInternal.peff  = flux_peff * delta_t;
    fluxInternal.euztw = flux_euztw * delta_t;
    fluxInternal.euzfw = flux_euzfw * delta_t;
    fluxInternal.elztw = flux_elztw * delta_t;
    fluxInternal.twexu = flux_twexu * delta_t;
    fluxInternal.ru    = flux_ru * delta_t;
    fluxInternal.pc    = flux_pc * delta_t;
    fluxInternal.pctw  = flux_pctw * delta_t;
    fluxInternal.pcfwp = flux_pcfwp * delta_t;
    fluxInternal.pcfws = flux_pcfws * delta_t;
    fluxInternal.twexl = flux_twexl * delta_t;
    fluxInternal.twexlp= flux_twexlp * delta_t;
    fluxInternal.twexls= flux_twexls * delta_t;
    fluxInternal.rlp   = flux_rlp * delta_t;
    fluxInternal.rls   = flux_rls * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




