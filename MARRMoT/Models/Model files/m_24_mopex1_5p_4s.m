function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_24_mopex1_5p_4s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: MOPEX-1
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Ye, S., Yaeger, M., Coopersmith, E., Cheng, L., & Sivapalan, M. (2012). 
% Exploring the physical controls of regional patterns of flow duration 
% curves - Part 2: Role of seasonality, the regime curve, and associated 
% process controls. Hydrology and Earth System Sciences, 16(11), 4447–4465.
% http://doi.org/10.5194/hess-16-4447-2012


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
tw      = theta(2);     % Groundwater leakage time [d-1]
tu      = theta(3);     % Slow flow routing response time [d-1]
se      = theta(4);     % Root zone storage capacity [mm]
tc      = theta(5);     % Mean residence time [d-1]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial soil moisture storage
S20     = storeInitial(2);       % Initial groundwater storage
S30     = storeInitial(3);       % Initial fast flow routing
S40     = storeInitial(4);       % Initial slow flow routing

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0];           % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);

flux_et1  = zeros(1,t_end);
flux_q1f  = zeros(1,t_end);
flux_qw   = zeros(1,t_end);
flux_et2  = zeros(1,t_end);
flux_q2u  = zeros(1,t_end);
flux_qf   = zeros(1,t_end);
flux_qs   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture
% S2. Groundwater
% S3. Fast flow routing
% S4. Slow flow routing

% ET1(S1,s1max,Ep(t),delta_t): evaporation from soil moisture
ET1 = evap_7;

% Q1F(P(t),S1,s1max): saturation excess
Q1F = saturation_1;

% QW(tw,S1): recharge to deeper store
QW = recharge_3;

% ET2(S2,se,Ep(t),delta_t): evaporation from deeper store
ET2 = evap_7;

% Q2U(tu,S2): slow response flow
Q2U = baseflow_1;

% QF(tc,S3): fast response routing
QF = baseflow_1;

% QS(tc,S4): fast response routing
QS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0;
                                               1,1,0,0;
                                               1,0,1,0;
                                               0,1,0,1]);                   % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0;
                                                  1,1,0,0;
                                                  1,0,1,0;
                                                  0,1,0,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    if t == 1; S4old = S40; else; S4old = store_S4(t-1); end                % store 4 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4) (P(t) - ...
                              ET1(S1,s1max,Ep(t),delta_t) - ...
                              Q1F(P(t),S1,s1max) - ...
                              QW(tw,S1));                                   % store 1 function with current flux values
    
    tmpf_S2 = @(S1,S2,S3,S4) (QW(tw,S1) - ...
                              ET2(S2,se,Ep(t),delta_t) - ...
                              Q2U(tu,S2));                                  % store 2 function
    
    tmpf_S3 = @(S1,S2,S3,S4) (Q1F(P(t),S1,s1max) - ...
                              QF(tc,S3));                                   % store 3 function
    
    tmpf_S4 = @(S1,S2,S3,S4) (Q2U(tu,S2) - ...
                              QS(tc,S4));                                   % store 4 function  
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4);                     % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4)),...        % system of storage equations
                        [S1old,S2old,S3old,S4old],...                       % storage values on previous time step
                        fsolve_options);                                    % solver options
    
    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       % [tmp_sNew,tmp_resnorm,flag]
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                                    eq_sys(1),eq_sys(2),...
                                                    eq_sys(3),eq_sys(4)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old], ...      % storages at previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_et1(t)  = ET1(tmp_sFlux(1),s1max,Ep(t),delta_t);
    flux_q1f(t)  = Q1F(P(t),tmp_sFlux(1),s1max);
    flux_qw(t)   = QW(tw,tmp_sFlux(1));
    flux_et2(t)  = ET2(tmp_sFlux(2),se,Ep(t),delta_t);
    flux_q2u(t)  = Q2U(tu,tmp_sFlux(2));
    flux_qf(t)   = QF(tc,tmp_sFlux(3));
    flux_qs(t)   = QS(tc,tmp_sFlux(4));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_et1(t) - flux_q1f(t) - flux_qw(t)) * delta_t;
    store_S2(t) = S2old + (flux_qw(t) - flux_et2(t) - flux_q2u(t)) * delta_t;    
    store_S3(t) = S3old + (flux_q1f(t) - flux_qf(t)) * delta_t;    
    store_S4(t) = S4old + (flux_q2u(t) - flux_qs(t)) * delta_t;    
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_et1 + flux_et2) * delta_t;
    fluxOutput.Q      = (flux_qf + flux_qs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.et1  = flux_et1 * delta_t;
    fluxInternal.et2  = flux_et2 * delta_t;
    fluxInternal.q1f  = flux_q1f * delta_t;
    fluxInternal.qw   = flux_qw * delta_t;
    fluxInternal.q2u  = flux_q2u * delta_t;
    fluxInternal.qf   = flux_qf * delta_t;
    fluxInternal.qs   = flux_qs * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

                       % delta S


 




