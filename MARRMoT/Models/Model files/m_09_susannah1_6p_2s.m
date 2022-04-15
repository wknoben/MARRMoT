function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_09_susannah1_6p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Susannah Brook model v1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Son, K., & Sivapalan, M. (2007). Improving model structure and reducing 
% parameter uncertainty in conceptual water balance models through the use 
% of auxiliary data. Water Resources Research, 43(1). 
% https://doi.org/10.1029/2006WR005032

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
sb  = theta(1);     % Maximum soil moisture storage [mm]
sfc = theta(2);     % Wiliting point as fraction of sb [-]
m   = theta(3);     % Fraction forest [-]
a   = theta(4);     % Runoff coefficient [d]  
b   = theta(5);     % Runoff coefficient [-]
r   = theta(6);     % Runoff coefficient [d-1]

%%INITIALISE MODEL STORES 
S10     = storeInitial(1);       % Initial upper zone storage
S20     = storeInitial(2);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                   % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_ebs  = zeros(1,t_end);
flux_eveg = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qss  = zeros(1,t_end);
flux_qr   = zeros(1,t_end);
flux_qb   = zeros(1,t_end);
flux_qt   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Upper zone
% S2. Groundwater

% EBS(m,S1,sb,Ep(t),delta_t): bare soil evaporation from soil moisture. Angle discontinuity
EBS = evap_5;

% EVEG(m,sfc,S1,sb,Ep(t),delta_t): vegetation transpiration
EVEG = evap_6;

% QSE(P,S1,sb): saturation excess overland flow
QSE = saturation_1;

% QSS(S1,sb,sfc,a,b,delta_t): non-linear sub-surface flow
QSS = interflow_7;

% QR(r,QSS): groundwater recharge
QR = baseflow_1;

% QB(S2,a,b,delta_t): baseflow
QB = baseflow_2;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % Disable display settings
                                'JacobPattern', [1,0;
                                                 1,1]);                     % Specify the Jacobian pattern                                             
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
    tmpf_S1 = @(S1,S2) (P(t) - ...
                        EBS(m,S1,sb,Ep(t),delta_t) - ...
                        EVEG(m,sfc,S1,sb,Ep(t),delta_t) - ...
                        QSE(P(t),S1,sb) - ...
                        QSS(S1,sb,sfc,a,b,delta_t));                        % Store 1 function

    tmpf_S2 = @(S1,S2) (QR(r,QSS(S1,sb,sfc,a,b,delta_t)) - ...
                        QB(S2,a,b,delta_t));                                % Store 2 function
    
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
    flux_ebs(t)  = EBS(m,tmp_sFlux(1),sb,Ep(t),delta_t);
    flux_eveg(t) = EVEG(m,sfc,tmp_sFlux(1),sb,Ep(t),delta_t);
    flux_qse(t)  = QSE(P(t),tmp_sFlux(1),sb);
    flux_qss(t)  = QSS(tmp_sFlux(1),sb,sfc,a,b,delta_t);
    flux_qr(t)   = QR(r,flux_qss(t));
    flux_qb(t)   = QB(tmp_sFlux(2),a,b,delta_t);
    flux_qt(t)   = flux_qse(t) + (flux_qss(t)-flux_qr(t)) + flux_qb(t);
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ebs(t) - flux_eveg(t) - flux_qse(t) - flux_qss(t)) * delta_t;
    store_S2(t) = S2old + (flux_qr(t) -  flux_qb(t)) * delta_t;

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = (flux_ebs + flux_eveg) * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.qr   = flux_qr * delta_t;
    fluxInternal.Ebs  = flux_ebs * delta_t;
    fluxInternal.Evg  = flux_eveg * delta_t;
    fluxInternal.Qse  = flux_qse * delta_t;
    fluxInternal.Qss  = flux_qss * delta_t;
    fluxInternal.Qb   = flux_qb * delta_t;    

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end



 




