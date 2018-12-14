function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_06_alpine1_4p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Alpine model v1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Eder, G., Sivapalan, M., & Nachtnebel, H. P. (2003). Modelling water 
% balances in an Alpine catchment through exploitation of emergent 
% properties over changing time scales. Hydrological Processes, 17(11), 
% 2125–2149. http://doi.org/10.1002/hyp.1325

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
tt   = theta(1);     % Threshold temperature for snowfall/snowmelt [celsius]
ddf  = theta(2);     % Degree-day-factor for snowmelt [mm/d/celsius]
Smax = theta(3);     % Maximum soil moisture storage [mm]
tc   = theta(4);     % Runoff coefficient [d-1]

%%INITIALISE MODEL STORES
S10  = storeInitial(1);       % Initial snow storage
S20  = storeInitial(2);       % Initial soil moisture storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];            % lower bounds of stores
store_upp = [];               % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_ps  = zeros(1,t_end);
flux_pr  = zeros(1,t_end);
flux_qn  = zeros(1,t_end);
flux_ea  = zeros(1,t_end);
flux_qse = zeros(1,t_end);
flux_qss = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Snow 
% S2. Soil moisture

% PS(P(t),T(t),tt): precipitation as snowfall
PS = snowfall_1;

% PR(P(t),T(t),tt): precipitation as rainfall
PR = rainfall_1;

% QN(ddf,tt,T(t),S1,delta_t): snowmelt
QN = melt_1;

% EA(S2,Ep(t),delta_t): evaporation from soil moisture
EA = evap_1;

% QSE(P(t)+QN(ddf,tt,T(t),S1),S2,Smax): saturation excess runoff
QSE = saturation_1;

% QSS(tc,S2): subsurface runoff
QSS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0;
                                               1,1]);                       % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0;
                                                  1,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                 % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                 % store 2 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2) (PS(P(t),T(t),tt) - ...
                        QN(ddf,tt,T(t),S1,delta_t));                     % store 1 function with current flux values
    tmpf_S2 = @(S1,S2) (PR(P(t),T(t),tt) + ...
                        QN(ddf,tt,T(t),S1,delta_t) - ...
                        EA(S2,Ep(t),delta_t) - ...
                        QSE(P(t)+QN(ddf,tt,T(t),S1,delta_t),S2,Smax) - ...
                        QSS(tc,S2));                                        % store 2 function
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2)),...                            % system of storage equations
                        [S1old,S2old],...                                   % storage values on previous time step
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
                                                    eq_sys(1),eq_sys(2)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old], ...                  % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ps(t)  = PS(P(t),T(t),tt);
    flux_pr(t)  = PR(P(t),T(t),tt);
    flux_qn(t)  = QN(ddf,tt,T(t),tmp_sFlux(1),delta_t);
    flux_ea(t)  = EA(tmp_sFlux(2),Ep(t),delta_t);
    flux_qse(t) = QSE(P(t)+flux_qn(t),tmp_sFlux(2),Smax);
    flux_qss(t) = QSS(tc,tmp_sFlux(2));
    
    % Update the stores
    store_S1(t) = S1old + (flux_ps(t) - flux_qn(t)) * delta_t;
    store_S2(t) = S2old + (flux_pr(t) + flux_qn(t) - flux_ea(t) - flux_qse(t) - flux_qss(t)) * delta_t;    

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = (flux_qse + flux_qss) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ps   = flux_ps * delta_t;
    fluxInternal.pr   = flux_pr * delta_t;
    fluxInternal.qn   = flux_qn * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qss  = flux_qss * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end





