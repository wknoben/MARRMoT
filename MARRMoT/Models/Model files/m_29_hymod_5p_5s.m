function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
    m_29_hymod_5p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: HyMOD
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Boyle, D. P. (2001). Multicriteria calibration of hydrologic models. 
% University of Arizona. Retrieved from http://hdl.handle.net/10150/290657
%
% Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, Hoshin, 
% V., & Sorooshian, S. (2001). A framework for development and application 
% of hydrological models. Hydrology and Earth System Sciences, 5, 13–26.

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
smax   = theta(1);     % Maximum soil moisture storage     [mm], 
b      = theta(2);     % Soil depth distribution parameter [-]
a      = theta(3);     % Runoff distribution fraction [-]
kf     = theta(4);     % Fast runoff coefficient [d-1]
ks     = theta(5);     % Slow runoff coefficient [d-1]

%%INITIALISE MODEL STORES AND ROUTING VECTORS
S10     = storeInitial(1);       % Initial soil moisture storage
S20     = storeInitial(2);       % Initial fast flow 1
S30     = storeInitial(3);       % Initial fast flow 2
S40     = storeInitial(4);       % Initial fast flow 3
S50     = storeInitial(5);       % Initial slow flow

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0];         % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS 
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);

flux_ea  = zeros(1,t_end);
flux_pe  = zeros(1,t_end);
flux_pf  = zeros(1,t_end);
flux_ps  = zeros(1,t_end);
flux_qf1 = zeros(1,t_end);
flux_qf2 = zeros(1,t_end);
flux_qf3 = zeros(1,t_end);
flux_qs  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture
% S2. Fast 1
% S3. Fast 2
% S4. Fast 3
% S5. Slow

% EA(S1,smax,Ep(t),delta_t): evap from soil moisture
EA = evap_7;

% PE(S1,smax,b,P(t))
PE = saturation_2;

% PF(a,PE(S1,smax,b,P(t))): split to fast reservoirs
PF = split_1;

% PS(1-a,PE(S1,smax,b,P(t))): split to slow reservoir
PS = split_1;

% QF1(kf,S2)
QF1 = baseflow_1;

% QF2(kf,S3)
QF2 = baseflow_1;

% QF3(kf,S4)
QF3 = baseflow_1;

% QS(ks,S5)
QS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...             % Disable display settings
                                'JacobPattern', [1,0,0,0,0;
                                                 1,1,0,0,0;
                                                 0,1,1,0,0;
                                                 0,0,1,1,0;
                                                 1,0,0,0,1]);           % Specify the Jacobian pattern                                             
lsqnonlin_options = optimoptions('lsqnonlin',...                        % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern',[1,0,0,0,0;
                                                 1,1,0,0,0;
                                                 0,1,1,0,0;
                                                 0,0,1,1,0;
                                                 1,0,0,0,1],...
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

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5) (P(t) - ...
                                 EA(S1,smax,Ep(t),delta_t) - ... 
                                 PE(S1,smax,b,P(t)));                       % Store 1 function with current flux values
    tmpf_S2 = @(S1,S2,S3,S4,S5) (PF(a,PE(S1,smax,b,P(t))) - ...
                                 QF1(kf,S2));                               % Store 2 funciton
    tmpf_S3 = @(S1,S2,S3,S4,S5) (QF1(kf,S2) - ...
                                 QF2(kf,S3));                               % Store 3 funciton
    tmpf_S4 = @(S1,S2,S3,S4,S5) (QF2(kf,S3) - ...
                                 QF3(kf,S4));                               % Store 4 funciton
    tmpf_S5 = @(S1,S2,S3,S4,S5) (PS(1-a,PE(S1,smax,b,P(t))) - ...
                                 QS(ks,S5));                                % Store 5 funciton
                             
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5);             % this returns a new anonymous function that we solve in the next step
                  
% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                        eq_sys(5)),...                                      % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old],...                 % storage values on previous time step
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
                                                    eq_sys(3),eq_sys(4),...
                                                    eq_sys(5)), ...
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
    flux_ea(t)  = EA(tmp_sFlux(1),smax,Ep(t),delta_t);
    flux_pe(t)  = PE(tmp_sFlux(1),smax,b,P(t));
    flux_pf(t)  = PF(a,flux_pe(t));
    flux_ps(t)  = PS(1-a,flux_pe(t));
    flux_qf1(t) = QF1(kf,tmp_sFlux(2));
    flux_qf2(t) = QF2(kf,tmp_sFlux(3));
    flux_qf3(t) = QF3(kf,tmp_sFlux(4));
    flux_qs(t)  = QS(ks,tmp_sFlux(5));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ea(t) - flux_pe(t)) * delta_t;
    store_S2(t) = S2old + (flux_pf(t) - flux_qf1(t)) * delta_t;
    store_S3(t) = S3old + (flux_qf1(t) - flux_qf2(t)) * delta_t;
    store_S4(t) = S4old + (flux_qf2(t) - flux_qf3(t)) * delta_t;
    store_S5(t) = S5old + (flux_ps(t) - flux_qs(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = (flux_qf3 + flux_qs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pe  = flux_pe * delta_t;
    fluxInternal.pf  = flux_pf * delta_t;
    fluxInternal.ps  = flux_ps * delta_t;
    fluxInternal.qf1 = flux_qf1 * delta_t;
    fluxInternal.qf2 = flux_qf2 * delta_t;
    fluxInternal.qf3 = flux_qf3 * delta_t;
    fluxInternal.qs  = flux_qs * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




