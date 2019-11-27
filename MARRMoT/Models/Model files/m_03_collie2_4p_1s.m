function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_03_collie2_4p_1s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Collie River Basin v2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Jothityangkoon, C., M. Sivapalan, and D. Farmer (2001), Process controls
% of water balance variability in a large semi-arid catchment: downward 
% approach to hydrological model development. Journal of Hydrology, 254,
% 174198. doi: 10.1016/S0022-1694(01)00496-6.


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
S1max   = theta(1);     % Maximum soil moisture storage [mm]
Sfc     = theta(2);     % Field capacity as fraction of S1max [-] 
a       = theta(3);     % Subsurface runoff coefficient [d-1]
M       = theta(4);     % Fraction forest cover [-]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial soil moisture storage

%%DEFINE STORE BOUNDARIES
store_min = [0];                     % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1  = zeros(1,t_end);
flux_eb   = zeros(1,t_end);
flux_ev   = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qss  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture

% EB(S1,S1max,(1-M)*Ep(t),delta_t): bare soil evaporation
EB = evap_7;

% EV(Sfc,S1,S1max,M*Ep(t),delta_t): evaporation from vegetation
EV = evap_3;

% QSE(P(t),S1,S1max): saturation excess overland flow
QSE = saturation_1;

% QSS(S1,a,Sfc*S1max): subsurface runoff
QSS = interflow_8;

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
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                 % store 1 at t-1
 
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1) (P(t) - ...
                     EB(S1,S1max,(1-M)*Ep(t),delta_t) - ...
                     EV(Sfc,S1,S1max,M*Ep(t),delta_t) - ...
                     QSE(P(t),S1,S1max) - ...
                     QSS(S1,a,Sfc*S1max));                                  % store 1 function with current flux values
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old],...
                      delta_t,...
                      tmpf_S1);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fzero(solve_fun,...
                                S1old,...
                                fzero_options);                             % 1 store solver

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
    flux_eb(t)   = EB(tmp_sFlux(1),S1max,(1-M)*Ep(t),delta_t);
    flux_ev(t)   = EV(Sfc,tmp_sFlux(1),S1max,M*Ep(t),delta_t);
    flux_qse(t)  = QSE(P(t),tmp_sFlux(1),S1max);
    flux_qss(t)  = QSS(tmp_sFlux(1),a,Sfc*S1max);
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_eb(t) - flux_ev(t) - flux_qse(t) - flux_qss(t)) * delta_t; 
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eb + flux_ev) * delta_t;
    fluxOutput.Q      = (flux_qse + flux_qss)  * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.eb   = flux_eb * delta_t;
    fluxInternal.ev   = flux_ev * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qss  = flux_qss * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end





