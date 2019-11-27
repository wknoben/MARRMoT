function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_21_flexb_9p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Flex-B
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
s1max   = theta(1);     % Maximum soil moisture storage [mm]
beta    = theta(2);     % Unsaturated zone shape parameter [-]
d       = theta(3);     % Fast/slow runoff distribution parameter [-]
percmax = theta(4);     % Maximum percolation rate [mm/d]
lp      = theta(5);     % Wilting point as fraction of s1max [-]
nlagf   = theta(6);     % Flow delay before fast runoff [d]
nlags   = theta(7);     % Flow delay before slow runoff [d]
kf      = theta(8);     % Fast runoff coefficient [d-1]
ks      = theta(9);     % Slow runoff coefficient [d-1]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial soil moisture storage
S20     = storeInitial(2);       % Initial fast routing storage
S30     = storeInitial(3);       % Initial slow routing storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];                 % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);

flux_ru  = zeros(1,t_end);
flux_eur = zeros(1,t_end);
flux_ps  = zeros(1,t_end);
flux_rf  = zeros(1,t_end);
flux_rs  = zeros(1,t_end);
flux_rfl = zeros(1,t_end);
flux_rsl = zeros(1,t_end);
flux_qf  = zeros(1,t_end);
flux_qs  = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
% [Optional]
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
% S1. Soil moisture
% S2. Fast reservoir
% S3. Slow reservoir

% RU(S1,s1max,beta,P(t)): infiltration into soil moisture
RU = saturation_3;

% EUR(lp,S1,s1max,Ep(t),delta_t): evaporation with wilting point
EUR = evap_3;

% PS(percmax,S1,s1max,delta_t): percolation to slow reservoir
PS = percolation_2;

% RS(d,P(t)-RU(S1,s1max,beta,P(t))): fraction non-infiltrated to slow flow
RS = split_1;

% RF(1-d,P(t)-RU(S1,s1max,beta,P(t))): fraction non-infiltrated to fast flow
RF = split_1;

% QF(kf,S2): fast flow
QF = baseflow_1;

% QS(ks,S3): fast flow
QS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
%   Store 1
fzero_options = optimset('Display','off');                                  % [1 store] settings of the root finding method
lsqnonlin_options1 = optimoptions('lsqnonlin',...                           % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'MaxFunEvals',1000);

%   Store 2 and 3
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0;
                                               0,1]);                       % Specify the Jacobian pattern                                               
lsqnonlin_options23 = optimoptions('lsqnonlin',...                          % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0;
                                                  0,1],...
                                 'MaxFunEvals',1000);
                             
%% 5. Solve the system for the full time series
for t = 1:t_end

% Special note: this model has two lag functions that delay the flows into 
% the fast and slow reservoirs. Due to the current setup of the routing
% functions, store 1 (soil moisture) needs to be computed first. The new S1
% value is then used for both routing schemes to determine inflow into
% store 2 (fast reservoir) and store 3(slow reservoir) respectively.
    
% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                 % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                 % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                 % store 2 at t-1
    
    % --- Solve store 1 (soil moisture) -----------------------------------
        % Create temporary store ODE's that need to be solved
        tmpf_S1 = @(S1) (RU(S1,s1max,beta,P(t)) - ...
                         EUR(lp,S1,s1max,Ep(t),delta_t) - ...
                         PS(percmax,S1,s1max,delta_t));
                     
        % Call the numerical scheme function to create the ODE approximations
        solve_fun1 = feval(scheme,...
                        [S1old],...
                        delta_t,...
                        tmpf_S1);        
        
        % Use the specified numerical scheme to solve storages
        [tmp_sNew(1), tmp_fval] = fzero(solve_fun1,...
                                    S1old,...
                                    fzero_options);                         % 1 store solver  
    
        % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(1),~,~] = rerunSolver('lsqnonlin', ...                
                                            lsqnonlin_options1, ...         % solver options
                                            @(eq_sys) solve_fun(...         % system of ODEs
                                                        eq_sys(1)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(1), ...                % recent estimates
                                            [S1old], ...                    % storages at previous time step
                                            store_min(1), ...               % lower bounds
                                            store_upp);                     % upper bounds              
        end
    
    % --- Handle the flow routing based on store 1 estimates --------------    
        % Find the storages needed to update fluxes: update 'tmp_sFlux'
        eval(store_fun);                                                    % creates/updates tmp_sFlux 

        % Fast routing
        tmp_rf_cur    = RF(1-d,P(t)-RU(tmp_sFlux(1),s1max,beta,P(t))).*uh_f;    % find how the runoff of this time step will be distributed in time
        tmp_rf_old    = tmp_rf_old + tmp_rf_cur;                                % update the 'still-to-flow-out' vector
        flux_rfl(t)   = tmp_rf_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_rf_old    = circshift(tmp_rf_old,-1);                               % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_rf_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)

        % --- Slow routing
        tmp_rs_cur    = (PS(percmax,tmp_sFlux(1),s1max,delta_t) + ...
                            RS(d,P(t)-RU(tmp_sFlux(1),s1max,beta,P(t)))).*uh_s; % find how the runoff of this time step will be distributed in time
        tmp_rs_old    = tmp_rs_old + tmp_rs_cur;                                % update the 'still-to-flow-out' vector
        flux_rsl(t)   = tmp_rs_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_rs_old    = circshift(tmp_rs_old,-1);                               % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_rs_old(end) = 0;    
    
    % --- Solve store 2 and 3 (fast and slow reservoirs) ------------------
        % Create temporary store ODE's that need to be solved
        tmpf_S2 = @(S2,S3) (flux_rfl(t) - ...
                            QF(kf,S2));
        tmpf_S3 = @(S2,S3) (flux_rsl(t) - ...
                            QS(ks,S3));
    
        % Call the numerical scheme function to create the ODE approximations
        solve_fun23 = feval(scheme,...
                          [S2old,S3old],...
                          delta_t,...
                          tmpf_S2,tmpf_S3);        

        % Use the specified numerical scheme to solve storages
        [tmp_sNew(2:3),tmp_fval] = fsolve(@(eq_sys) solve_fun23(...
                            eq_sys(1),eq_sys(2)),...                        % system of storage equations
                            [S2old,S3old],...                               % storage values on previous time step
                            fsolve_options);                                % solver options
 
        % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(2:3),~,~] = rerunSolver('lsqnonlin', ...              
                                            lsqnonlin_options23, ...        % solver options
                                            @(eq_sys) solve_fun23(...       % system of ODEs
                                            eq_sys(1),eq_sys(2)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(2:3), ...              % recent estimates
                                            [S2old,S3old], ...              % storages at previous time step
                                            store_min(2:3), ...             % lower bounds
                                            store_upp);                     % upper bounds              
        end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

   % Calculate the fluxes
    flux_ru(t)  = RU(tmp_sFlux(1),s1max,beta,P(t));
    flux_eur(t) = EUR(lp,tmp_sFlux(1),s1max,Ep(t),delta_t);
    flux_ps(t)  = PS(percmax,tmp_sFlux(1),s1max,delta_t);
    flux_rf(t)  = RF(1-d,P(t)-RU(tmp_sFlux(1),s1max,beta,P(t)));
    flux_rs(t)  = RS(d,P(t)-RU(tmp_sFlux(1),s1max,beta,P(t)));
    flux_qf(t)  = QF(kf,tmp_sFlux(2));
    flux_qs(t)  = QS(ks,tmp_sFlux(3));
    
    % Update the stores
    store_S1(t) = S1old + (flux_ru(t) - flux_eur(t) - flux_ps(t)) * delta_t;
    store_S2(t) = S2old + (flux_rfl(t) - flux_qf(t)) * delta_t;    
    store_S3(t) = S3old + (flux_rsl(t) - flux_qs(t)) * delta_t; 
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_eur * delta_t;
    fluxOutput.Q      = (flux_qf + flux_qs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ru   = flux_ru * delta_t;
    fluxInternal.ps   = flux_ps * delta_t;
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
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     tmp_rout.*delta_t);% Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end




