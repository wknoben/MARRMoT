function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
        m_04_newzealand1_6p_1s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: New Zealand model v1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Atkinson, S. E., Woods, R. A., & Sivapalan, M. (2002). Climate and
% landscape controls on water balance model complexity over changing 
% timescales. Water Resources Research, 38(12), 17–50. 
% http://doi.org/10.1029/2002WR001487


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
sfc     = theta(2);     % Field capacity as fraction of maximum soil moisture [-]
m       = theta(3);     % Fraction forest [-]
a       = theta(4);     % Subsurface runoff coefficient [d-1]
b       = theta(5);     % Runoff non-linearity [-]
tcbf    = theta(6);     % Baseflow runoff coefficient [d-1]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);          % Initial soil moisture storage

%%DEFINE STORE BOUNDARIES
store_min = [0];                    % lower bounds of stores
store_upp = [];                     % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);

flux_veg  = zeros(1,t_end);
flux_ebs  = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qss  = zeros(1,t_end);
flux_qbf  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture

% EVEG(m,sfc,S1,s1max,Ep(t),delta_t): transpiration through vegetation
EVEG = evap_6;

% EBS(m,S1,s1max,Ep(t),delta_t): bare soil evaporation
EBS = evap_5;

% QSE(P(t),S1,s1max): saturation excess overland flow
QSE = saturation_1;

% QSS(S1,a,sfc*s1max,b,delta_t): subsurface flow
QSS = interflow_9;

% QBF(tcbf,S1): baseflow. Angle discontinuity at S=0
QBF = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options                                                      
fzero_options = optimset('Display','off');                                  % [1 store] settings of the root finding method
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
   
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1) (P(t) - ...
                     EVEG(m,sfc,S1,s1max,Ep(t),delta_t) - ...
                     EBS(m,S1,s1max,Ep(t),delta_t) - ...
                     QSE(P(t),S1,s1max) - ...
                     QSS(S1,a,sfc*s1max,b,delta_t) - ...
                     QBF(tcbf,S1));                                         % store 1 function with current flux values
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old],...
                      delta_t,...
                      tmpf_S1);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew, tmp_fval] = fzero(solve_fun,...
                                    S1old,...
                                    fzero_options);                         % 1 store solver

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       % [tmp_sNew,tmp_resnorm,flag]
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                                    eq_sys(1)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old], ...                        % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_veg(t)  = EVEG(m,sfc,tmp_sFlux(1),s1max,Ep(t),delta_t);
    flux_ebs(t)  = EBS(m,tmp_sFlux(1),s1max,Ep(t),delta_t);
    flux_qse(t)  = QSE(P(t),tmp_sFlux(1),s1max);
    flux_qss(t)  = QSS(tmp_sFlux(1),a,sfc*s1max,b,delta_t);
    flux_qbf(t)  = QBF(tcbf,tmp_sFlux(1));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) -flux_veg(t) -flux_ebs(t) -flux_qse(t) -flux_qss(t) -flux_qbf(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_veg + flux_ebs) * delta_t;
    fluxOutput.Q      = (flux_qse + flux_qss + flux_qbf) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.eveg = flux_veg * delta_t;
    fluxInternal.ebs  = flux_ebs * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qss  = flux_qss * delta_t;
    fluxInternal.qbf  = flux_qbf * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    
% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end



 




