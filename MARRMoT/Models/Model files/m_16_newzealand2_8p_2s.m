function [ fluxOutput, fluxInternal, storeInternal, waterBalance  ] = ...
        m_16_newzealand2_8p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: New Zealand model v2 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Atkinson, S. E., Sivapalan, M., Woods, R. A., & Viney, N. R. (2003). 
% Dominant physical controls on hourly flow predictions and the role of 
% spatial variability: Mahurangi catchment, New Zealand. Advances in Water 
% Resources, 26(3), 219–235. http://doi.org/10.1016/S0309-1708(02)00183-5


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
s1max   = theta(1);     % Maximum interception storage [mm] 
s2max   = theta(2);     % Maximum soil moisture storage [mm] 
sfc     = theta(3);     % Field capacity as fraction of maximum soil moisture [-]
m       = theta(4);     % Fraction forest [-]
a       = theta(5);     % Subsurface runoff coefficient [d-1]
b       = theta(6);     % Runoff non-linearity [-]
tcbf    = theta(7);     % Baseflow runoff coefficient [d-1]
d       = theta(8);     % Routing time delay [d]

%%INITIALISE MODEL STORES
S10         = storeInitial(1);       % Initial interception storage
S20         = storeInitial(2);       % Initial soil moisture storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                   % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_eint = zeros(1,t_end);
flux_qtf  = zeros(1,t_end);
flux_veg  = zeros(1,t_end);
flux_ebs  = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qss  = zeros(1,t_end);
flux_qbf  = zeros(1,t_end);
flux_qt   = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPH
[~,uh_full] = uh_4_full(1,d,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_Qt_old  = zeros(1,length(uh_full));      % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Soil moisture

% EINT(S1,Ep(t),delta_t): evaporation from interception
EINT = evap_1;

% QTF(P(t),S1,s1max):
QTF = interception_1;

% EVEG(m,sfc,S2,s2max,Ep(t),delta_t): transpiration through vegetation
EVEG = evap_6;

% EBS(m,S2,s2max,Ep(t),delta_t): bare soil evaporation
EBS = evap_5;

% QSE(QTF(P(t),S1,s1max),S2,s2max): saturation excess overland flow
QSE = saturation_1;

% QSS(S2,a,sfc*s2max,b,delta_t): subsurface flow
QSS = interflow_9;

% QBF(tcbf,S2): baseflow
QBF = baseflow_1;

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
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2) (P(t) - ...
                        EINT(S1,Ep(t),delta_t) - ...
                        QTF(P(t),S1,s1max));                                % store 1 function with current flux values
    tmpf_S2 = @(S1,S2) (QTF(P(t),S1,s1max) - ...
                        EVEG(m,sfc,S2,s2max,Ep(t),delta_t) - ...
                        EBS(m,S2,s2max,Ep(t),delta_t) - ...
                        QSE(QTF(P(t),S1,s1max),S2,s2max) - ...
                        QSS(S2,a,sfc*s2max,b,delta_t) - ...
                        QBF(tcbf,S2));                                      % store 2 function
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2);                                     % this returns a new anonymous function that we solve in the next step

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
    % Calculate the fluxes
    flux_eint(t) = EINT(tmp_sNew(1),Ep(t),delta_t);
    flux_qtf(t)  = QTF(P(t),tmp_sNew(1),s1max);
    flux_veg(t)  = EVEG(m,sfc,tmp_sNew(2),s2max,Ep(t),delta_t);
    flux_ebs(t)  = EBS(m,tmp_sNew(2),s2max,Ep(t),delta_t);
    flux_qse(t)  = QSE(QTF(P(t),tmp_sNew(1),s1max),tmp_sNew(2),s2max);
    flux_qss(t)  = QSS(tmp_sNew(2),a,sfc*s2max,b,delta_t);
    flux_qbf(t)  = QBF(tcbf,tmp_sNew(2));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_eint(t) - flux_qtf(t)) * delta_t;
    store_S2(t) = S2old + (flux_qtf(t) -flux_veg(t) -flux_ebs(t) - ...
                           flux_qse(t) -flux_qss(t) -flux_qbf(t)) * delta_t;

    % Apply the routing scheme
    tmp_Qt_cur    = (flux_qse(t) + flux_qss(t) + flux_qbf(t)).*uh_full;     % find how the runoff of this time step will be distributed in time
    tmp_Qt_old    = tmp_Qt_old + tmp_Qt_cur;                                % update the 'still-to-flow-out' vector
    flux_qt(t)    = tmp_Qt_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_Qt_old = circshift(tmp_Qt_old,-1);                                  % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_Qt_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eint + flux_veg + flux_ebs) * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.eint = flux_eint * delta_t;
    fluxInternal.qtf  = flux_qtf * delta_t;
    fluxInternal.eveg = flux_veg * delta_t;
    fluxInternal.ebs  = flux_ebs * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qss  = flux_qss * delta_t;
    fluxInternal.qbf  = flux_qbf * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     tmp_Qt_old);       % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


 




