function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_34_flexis_12p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: FLEX-IS 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Fenicia, F., McDonnell, J. J., & Savenije, H. H. G. (2008). Learning from
% model improvement: On the contribution of complementary data to process
% understanding. Water Resources Research, 44(6), 1–13. 
% http://doi.org/10.1029/2007WR006386
%
% Nijzink, R., Hutton, C., Pechlivanidis, I., Capell, R., Arheimer, B., 
% Freer, J., … Hrachowitz, M. (2016). The evolution of root zone moisture 
% capacities after land use change: a step towards predictions under 
% change? Hydrology and Earth System Sciences Discussions, 20(August), 
% 4775–4799. http://doi.org/10.5194/hess-2016-427


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

% Parameters 
% [name in documentation] = theta(order in which specified in parameter file)
smax    = theta(1);     % Maximum soil moisture storage [mm]
beta    = theta(2);     % Unsaturated zone shape parameter [-]
d       = theta(3);     % Fast/slow runoff distribution parameter [-]
percmax = theta(4);     % Maximum percolation rate [mm/d]
lp      = theta(5);     % Wilting point as fraction of s1max [-]
nlagf   = theta(6);     % Flow delay before fast runoff [d]
nlags   = theta(7);     % Flow delay before slow runoff [d]
kf      = theta(8);     % Fast runoff coefficient [d-1]
ks      = theta(9);     % Slow runoff coefficient [d-1]
imax    = theta(10);    % Maximum interception storage [mm]
tt      = theta(11);    % Threshold temperature for snowfall/snowmelt [oC]
ddf     = theta(12);    % Degree-day factor for snowmelt [mm/d/oC]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial snow pack
S20     = storeInitial(2);       % Initial interception
S30     = storeInitial(3);       % Initial soil moisture storage
S40     = storeInitial(4);       % Initial fast store
S50     = storeInitial(5);       % Initial slow store

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0];             % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);

flux_ps  = zeros(1,t_end);
flux_pi  = zeros(1,t_end);
flux_m   = zeros(1,t_end);
flux_peff= zeros(1,t_end);
flux_ei  = zeros(1,t_end);
flux_ru  = zeros(1,t_end);
flux_eur = zeros(1,t_end);
flux_rp  = zeros(1,t_end);
flux_rf  = zeros(1,t_end);
flux_rs  = zeros(1,t_end);
flux_rfl = zeros(1,t_end);
flux_rsl = zeros(1,t_end);
flux_qf  = zeros(1,t_end);
flux_qs  = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_f] = uh_3_half(1,nlagf,delta_t);
[~,uh_s] = uh_3_half(1,nlags,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_rf_old  = zeros(1,length(uh_f));  % temporary vector needed to deal with routing
tmp_rs_old  = zeros(1,length(uh_s));  % temporary vector needed to deal with routing

%%INITIALISE FLUX UPDATE VECTOR
tmp_sNew = NaN.*zeros(1,length(storeInitial));

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Snow pack
% S2. Interception
% S3. Soil moisture
% S4. Fast reservoir
% S5. Slow reservoir

% PS(P(t),T(t),tt): snowfall
PS = snowfall_1; 

% M(ddf,tt,T(t),S1,delta_t): snowmelt
M = melt_1;

% PI(P(t),T(t),tt): precipitation as rain
PI = rainfall_1;

% PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax): throughfall from 
% interception
PEFF = interception_1;

% EI(S2,Ep(t),delta_t): evaporation from interception
EI = evap_1;

% RU(S3,smax,beta,PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax)): 
% infiltration into soil moisture
RU = saturation_3;

% EUR(lp,S3,smax,Ep(t),delta_t): evaporation with wilting point
EUR = evap_3;

% RP(percmax,S3,smax,delta_t): percolation to slow reservoir
RP = percolation_2;

% RS(d,PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax)-
%  RU(S3,smax,beta,PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax))):
% fraction non-infiltrated to slow flow
RS = split_1;

% RF(1-d,PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax)-
%  RU(S3,smax,beta,PEFF( M(ddf,tt,T(t),S1,delta_t)+PP(P(t),T(t),tt) ,S2,imax))):
% fraction non-infiltrated to slow flow
RF = split_1;

% QF(kf,S4): fast flow
QF = baseflow_1;

% QS(ks,S5): fast flow
QS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
%   Store 1, 2 & 3 (snow, interception, soil moisture)
fsolve_options123 = optimoptions('fsolve','Display','none',...              % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0;
                                               1,1,0;
                                               1,1,1]);                     % Specify the Jacobian pattern                                               
lsqnonlin_options123 = optimoptions('lsqnonlin',...                         % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0;
                                                  1,1,0;
                                                  1,1,1],...
                                 'MaxFunEvals',1000);

%   Store 4 & 5 (fast and slow stores)
fsolve_options45 = optimoptions('fsolve','Display','none',...               % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0;
                                               0,1]);                       % Specify the Jacobian pattern                                               
lsqnonlin_options45 = optimoptions('lsqnonlin',...                          % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0;
                                                  0,1],...
                                 'MaxFunEvals',1000);
                             
%% 5. Solve the system for the full time series
for t = 1:t_end

% Special note: this model has two lag functions that delay the flows into 
% the fast and slow reservoirs. Due to the current setup of the routing
% functions, store 1 & 2 (interception and soil moisture) need to be 
% computed first. The new S2 value is then used for both routing schemes to 
% determine inflow into store 3 (fast reservoir) and store 4 (slow
% reservoir) respectively.    
    
% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    if t == 1; S4old = S40; else; S4old = store_S4(t-1); end                % store 4 at t-1
    if t == 1; S5old = S50; else; S5old = store_S5(t-1); end                % store 5 at t-1

    % --- Solve store 1, 2 AND 3 (interception and soil moisture) ---------
        % Create temporary store ODE's that need to be solved
        tmpf_S1 = @(S1,S2,S3) (PS(P(t),T(t),tt) - ...
                                M(ddf,tt,T(t),S1,delta_t));
        
        tmpf_S2 = @(S1,S2,S3) (PI(P(t),T(t),tt) + ...
                                M(ddf,tt,T(t),S1,delta_t) - ...
                                EI(S2,Ep(t),delta_t) - ...
                                PEFF( M(ddf,tt,T(t),S1,delta_t)+PI(P(t),T(t),tt) ,S2,imax));
        
        tmpf_S3 = @(S1,S2,S3) (RU(S3,smax,beta,PEFF( M(ddf,tt,T(t),S1,delta_t)+PI(P(t),T(t),tt) ,S2,imax)) - ...
                                EUR(lp,S3,smax,Ep(t),delta_t) - ...
                                RP(percmax,S3,smax,delta_t));
        
        % Call the numerical scheme function to create the ODE approximations
        solve_fun123 = feval(scheme,...
                        [S1old,S2old,S3old],...
                        delta_t,...
                        tmpf_S1,tmpf_S2,tmpf_S3);        
        
        % Use the specified numerical scheme to solve storages
        [tmp_sNew(1:3),tmp_fval] = fsolve(@(eq_sys) solve_fun123(...
                            eq_sys(1),eq_sys(2),eq_sys(3)),...              % system of storage equations
                            [S1old,S2old,S3old],...                         % storage values on previous time step
                            fsolve_options123);  
    
        % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(1:3),~,~] = rerunSolver('lsqnonlin', ...              
                                            lsqnonlin_options123, ...       % solver options
                                            @(eq_sys) solve_fun123(...      % system of ODEs
                                            eq_sys(1),eq_sys(2),eq_sys(3)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(1:3), ...              % recent estimates
                                            [S1old,S2old,S3old], ...        % storages ate previous time step
                                            store_min(1:3), ...             % lower bounds
                                            store_upp);                     % upper bounds              
        end
    
    % --- Handle the flow routing based on store 2 estimates --------------    
        % Find the storages needed to update fluxes: update 'tmp_sFlux'
        eval(store_fun);                                                    % creates/updates tmp_sFlux 

        % Fast routing
        tmp_rf_cur    = RF(1-d,PEFF( M(ddf,tt,T(t),tmp_sFlux(1),delta_t)+...
                        PI(P(t),T(t),tt) ,tmp_sFlux(2),imax)-...
                        RU(tmp_sFlux(3),smax,beta,PEFF( M(ddf,tt,T(t),tmp_sFlux(1),delta_t)+...
                        PI(P(t),T(t),tt) ,tmp_sFlux(2),imax))).*uh_f;       % find how the runoff of this time step will be distributed in time
        tmp_rf_old    = tmp_rf_old + tmp_rf_cur;                            % update the 'still-to-flow-out' vector
        flux_rfl(t)   = tmp_rf_old(1);                                      % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_rf_old    = circshift(tmp_rf_old,-1);                           % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_rf_old(end) = 0;                                                % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)

        % --- Slow routing
        tmp_rs_cur    = (RS(d,PEFF( M(ddf,tt,T(t),tmp_sFlux(1),delta_t)+...
                        PI(P(t),T(t),tt) ,tmp_sFlux(2),imax)-...
                        RU(tmp_sFlux(3),smax,beta,PEFF( M(ddf,tt,T(t),tmp_sFlux(1),delta_t)+...
                        PI(P(t),T(t),tt) ,tmp_sFlux(2),imax)))+...
                        RP(percmax,tmp_sFlux(3),smax,delta_t)).*uh_s;       % find how the runoff of this time step will be distributed in time
        tmp_rs_old    = tmp_rs_old + tmp_rs_cur;                            % update the 'still-to-flow-out' vector
        flux_rsl(t)   = tmp_rs_old(1);                                      % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_rs_old    = circshift(tmp_rs_old,-1);                           % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_rs_old(end) = 0;
        
    % --- Solve store 4 and 5 (fast and slow reservoirs) ------------------
        % Create temporary store ODE's that need to be solved
        tmpf_S4 = @(S4,S5) (flux_rfl(t) - ...
                            QF(kf,S4));
        
        tmpf_S5 = @(S4,S5) (flux_rsl(t) - ...
                            QS(ks,S5));    
    
    % Call the numerical scheme function to create the ODE approximations
        solve_fun45 = feval(scheme,...
                          [S4old,S5old],...
                          delta_t,...
                          tmpf_S4,tmpf_S5);
    
    % Use the specified numerical scheme to solve storages
        [tmp_sNew(4:5),tmp_fval] = fsolve(@(eq_sys) solve_fun45(...
                            eq_sys(1),eq_sys(2)),...                        % system of storage equations
                            [S4old,S5old],...                               % storage values on previous time step
                            fsolve_options45);                              % solver options
 
    % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(4:5),~,~] = rerunSolver('lsqnonlin', ...              
                                            lsqnonlin_options45, ...        % solver options
                                            @(eq_sys) solve_fun45(...       % system of ODEs
                                            eq_sys(1),eq_sys(2)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(4:5), ...              % recent estimates
                                            [S4old,S5old], ...              % storages ate previous time step
                                            store_min(4:5), ...             % lower bounds
                                            store_upp);                     % upper bounds              
        end
                        
                        
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ps(t)  = PS(P(t),T(t),tt);
    flux_pi(t)  = PI(P(t),T(t),tt);
    flux_m(t)   = M(ddf,tt,T(t),tmp_sFlux(1),delta_t);
    flux_peff(t)= PEFF( flux_m(t)+flux_pi(t) ,tmp_sFlux(2),imax);
    flux_ei(t)  = EI(tmp_sFlux(2),Ep(t),delta_t);
    flux_ru(t)  = RU(tmp_sFlux(3),smax,beta,flux_peff(t));
    flux_eur(t) = EUR(lp,tmp_sFlux(3),smax,Ep(t),delta_t);
    flux_rp(t)  = RP(percmax,tmp_sFlux(3),smax,delta_t);
    flux_rf(t)  = RF(1-d,flux_peff(t)-flux_ru(t));
    flux_rs(t)  = RS(d,flux_peff(t)-flux_ru(t));
    flux_qf(t)  = QF(kf,tmp_sFlux(4));
    flux_qs(t)  = QS(ks,tmp_sFlux(5));
    
    % Update the stores
    store_S1(t) = S1old + (flux_ps(t) - flux_m(t)) * delta_t;
    store_S2(t) = S2old + (flux_m(t) + flux_pi(t) - flux_peff(t) - flux_ei(t)) * delta_t;
    store_S3(t) = S3old + (flux_ru(t) - flux_eur(t) - flux_rp(t)) * delta_t;
    store_S4(t) = S4old + (flux_rfl(t) - flux_qf(t)) * delta_t;
    store_S5(t) = S5old + (flux_rsl(t) - flux_qs(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ei + flux_eur) * delta_t;
    fluxOutput.Q      = (flux_qf + flux_qs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ps   = flux_ps * delta_t;
    fluxInternal.pi   = flux_pi * delta_t;
    fluxInternal.m    = flux_m * delta_t;
    fluxInternal.ei   = flux_ei * delta_t;
    fluxInternal.eur  = flux_eur * delta_t;
    fluxInternal.peff = flux_peff * delta_t;
    fluxInternal.ru   = flux_ru * delta_t;
    fluxInternal.rp   = flux_rp * delta_t;
    fluxInternal.rf   = flux_rf * delta_t;
    fluxInternal.rs   = flux_rs * delta_t;
    fluxInternal.rfl  = flux_rfl * delta_t;
    fluxInternal.rsl  = flux_rsl * delta_t;
    fluxInternal.qf   = flux_qf * delta_t;
    fluxInternal.qs   = flux_qs * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;

% Check water balance
if nargout == 4
    
    % Combine the routing vectors into one for water balance calculation
    if numel(tmp_rf_old) > numel(tmp_rs_old)
        tmp_rout = tmp_rf_old + [tmp_rs_old, zeros(1,numel(tmp_rf_old)-numel(tmp_rs_old))];
    elseif numel(tmp_rf_old) < numel(tmp_rs_old)
        tmp_rout = tmp_rs_old + [tmp_rf_old, zeros(1,numel(tmp_rs_old)-numel(tmp_rf_old))];
    else
        tmp_rout = tmp_rf_old + tmp_rs_old;
    end
    
    % Do the water balance check
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     tmp_rout);         % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end
 




