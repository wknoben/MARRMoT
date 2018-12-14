function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_25_tcm_6p_4s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Thames Catchment Model
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Moore, R. J., & Bell, V. A. (2001). Comparison of rainfall-runoff models 
% for flood forecasting. Part 1: Literature review of models. Bristol: 
% Environment Agency.

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
phi   = theta(1);     % Fraction preferential recharge [-]
rc    = theta(2);     % Maximum soil moisture depth [mm]
gam   = theta(3);     % Fraction of Ep reduction with depth [-]
k1    = theta(4);     % Runoff coefficient [d-1]
fa    = theta(5);     % Fraction of mean(P) that forms abstraction rate [mm/d]
k2    = theta(6);     % Runoff coefficient [mm-1 d-1]

% Auxiliary parameters
ca    = fa*mean(P);   % Abstraction rate [mm/d]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial upper soil moisture storage
S20   = storeInitial(2);       % Initial lower soil moisture deficit storage
S30   = storeInitial(3);       % Initial unsaturated zone
S40   = storeInitial(4);       % Initial saturated zone

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0];         % lower bounds of stores
store_upp = [];                % optional higher bounds

%%INITIALISE STORAGE VECTORS 
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);

flux_en     = zeros(1,t_end);
flux_ea     = zeros(1,t_end);
flux_et     = zeros(1,t_end);
flux_pn     = zeros(1,t_end);
flux_pby    = zeros(1,t_end);
flux_pin    = zeros(1,t_end);
flux_qex1   = zeros(1,t_end);
flux_qex2   = zeros(1,t_end);
flux_quz    = zeros(1,t_end);
flux_a      = zeros(1,t_end);
flux_q      = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Upper soil moisture
% S2. Lower soil moisture deficit
% S3. Unsaturated zone
% S4. Saturated zone

% PN(P(t),Ep(t)): effective rainfall
PN = effective_1;

% PIN(1-phi,PN(P(t),Ep(t))): infiltrated net precip
PIN = split_1;

% PBY(phi,PN(P(t),Ep(t))): preferential recharge
PBY = split_1;

% EA(S1,Ep(t),delta_t): evaporation from upper store
EA = evap_1;

% QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc): overflow from upper store
QEX1 = saturation_1;

% ET(gam,S2,S1,0.01,Ep(t),delta_t): reduced evaporation from lower store 
ET = evap_16;

% QEX2(QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc),S2,0.01): overflow from lower soil moisture
QEX2 = saturation_9;

% QUZ(k1,S3): outflow from unsaturated zone
QUZ = baseflow_1;

% A(ca): groundwater abstraction
A = abstraction_1;

% Q(k2,0,S4,delta_t): groundwater outflow
Q = baseflow_6;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0;
                                               1,1,0,0;
                                               1,1,1,0
                                               0,0,1,1]);                   % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,0,0,0;
                                                    1,1,0,0;
                                                    1,1,1,0
                                                    0,0,1,1],...
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
    tmpf_S1 = @(S1,S2,S3,S4) ...
               (PIN(1-phi,PN(P(t),Ep(t))) - ...
                EA(S1,Ep(t),delta_t) - ...
                QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc));                     % store 1 function with current flux values
            
    tmpf_S2 = @(S1,S2,S3,S4) ...                                            % store 2: DEFICIT STORE
               (ET(gam,S2,S1,0.01,Ep(t),delta_t) + ...
                QEX2(QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc),S2,0.01) - ...
                QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc));
            
    tmpf_S3 = @(S1,S2,S3,S4) ...                                            % store 3 function
               (PBY(phi,PN(P(t),Ep(t))) + ...
                QEX2(QEX1(PIN(1-phi,PN(P(t),Ep(t))),S1,rc),S2,0.01) - ...
                QUZ(k1,S3));

    tmpf_S4 = @(S1,S2,S3,S4) ...                                            % store 4 function
               (QUZ(k1,S3) - ...
                A(ca) - ...
                Q(k2,0,S4,delta_t));
        
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
                                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old], ...      % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_pn(t)     = PN(P(t),Ep(t));
    flux_en(t)     = P(t) - flux_pn(t);
    flux_pby(t)    = PBY(phi,flux_pn(t));
    flux_pin(t)    = PIN(1-phi,flux_pn(t));
    flux_ea(t)     = EA(tmp_sFlux(1),Ep(t),delta_t);
    flux_et(t)     = ET(gam,tmp_sFlux(2),tmp_sFlux(1),0.01,Ep(t),delta_t);
    flux_qex1(t)   = QEX1(flux_pin(t),tmp_sFlux(1),rc);
    flux_qex2(t)   = QEX2(flux_qex1(t),tmp_sFlux(2),0.01);
    flux_quz(t)    = QUZ(k1,tmp_sFlux(3));
    flux_a(t)      = A(ca);
    flux_q(t)      = Q(k2,0,tmp_sFlux(4),delta_t);
    
    % Update the stores
    store_S1(t) = S1old + (flux_pin(t) - flux_ea(t) - flux_qex1(t)) * delta_t;
    store_S2(t) = S2old + (flux_et(t)  + flux_qex2(t) - flux_qex1(t)) * delta_t;    
    store_S3(t) = S3old + (flux_qex2(t) + flux_pby(t) - flux_quz(t)) * delta_t;
    store_S4(t) = S4old + (flux_quz(t) - flux_a(t) - flux_q(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_en + flux_ea + flux_et) * delta_t;
    fluxOutput.Q      = flux_q * delta_t;
    fluxOutput.A      = flux_a * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pn   = flux_pn * delta_t;
    fluxInternal.en   = flux_en * delta_t;
    fluxInternal.pby  = flux_pby * delta_t;
    fluxInternal.pin  = flux_pin * delta_t;
    fluxInternal.ea   = flux_ea * delta_t;
    fluxInternal.et   = flux_et * delta_t;
    fluxInternal.qex1 = flux_qex1 * delta_t;
    fluxInternal.qex2 = flux_qex2 * delta_t;
    fluxInternal.quz  = flux_quz * delta_t;
    fluxInternal.a    = flux_a * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;

% Check water balance
% NOTE: this function doesn't yet work with deficit stores
% or abstractions. Thus manual calculation is needed
if nargout == 4
    waterBalance = ...
            sum(P) - ...                    % Incoming precipitation
            sum(fluxOutput.Ea) - ...        % Outgoing evaporation
            sum(fluxOutput.A) - ...         % Outgoing abstractions
            sum(fluxOutput.Q) - ...         % Outgoing runoff
            (store_S1(end)-S10) + ...       % Store change
            (store_S2(end)-S20) - ...       % Store deficit change
            (store_S3(end)-S30) - ...       % Store change
            (store_S4(end)-S40);            % Store change

    disp(['Total P  = ',num2str(sum(P)),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Total A  = ',num2str(sum(fluxOutput.A)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm deficit.'])
    disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
    disp(['Delta S4 = ',num2str((store_S4(end)-S40)),'mm.'])
    disp(['Water balance = sum(P) - sum(Ea) - sum(Q) - sum(A) - delta S = ',...
        num2str(waterBalance)]);                         
end


 




