function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_28_xinanjiang_12p_4s( fluxInput, storeInitial, theta, solver)
% Hydrologic conceptual model: Xinanjiang
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Jayawardena, A. W., & Zhou, M. C. (2000). A modified spatial soil moisture
% storage capacity distribution curve for the Xinanjiang model. Journal of 
% Hydrology, 227(1-4), 93–113. http://doi.org/10.1016/S0022-1694(99)00173-0
% 
% Zhao, R.-J. (1992). The Xinanjiang model applied in China. Journal of 
% Hydrology, 135(1-4), 371–381. http://doi.org/10.1016/0022-1694(92)90096-E


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
aim  = theta(1);     % Fraction impervious area [-]
a    = theta(2);     % Tension water distribution inflection parameter [-]
b    = theta(3);     % Tension water distribution shape parameter [-]
stot = theta(4);     % Total soil moisture storage (W+S) [mm]
fwmx = theta(5);     % Fraction of Stot that is Wmax [-]
flm  = theta(6);     % Fraction of wmax that is LM [-]
c    = theta(7);     % Fraction of LM for second evaporation change [-]
ex   = theta(8);     % Free water distribution shape parameter [-]
ki   = theta(9);     % Free water interflow parameter [d-1]
kg   = theta(10);    % Free water baseflow parameter [d-1]
ci   = theta(11);    % Interflow time coefficient [d-1]
cg   = theta(12);    % Baseflow time coefficient [d-1]

% Auxiliary parameters
wmax = fwmx*stot;    % Maximum tension water depth [mm]
smax = (1-fwmx)*stot;% Maximum free water depth [mm]
lm   = flm*wmax;     % Tension water threshold for evaporation change [mm]

%%INITIALISE MODEL STORES
S10  = storeInitial(1);       % Initial tension water storage
S20  = storeInitial(2);       % Initial free water storage
S30  = storeInitial(3);       % Initial interflow routing
S40  = storeInitial(4);       % Initial baseflow routing

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0];        % lower bounds of stores
store_upp = [];               % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);

flux_rb  = zeros(1,t_end);
flux_pi  = zeros(1,t_end);
flux_e   = zeros(1,t_end);
flux_r   = zeros(1,t_end);
flux_rs  = zeros(1,t_end);
flux_ri  = zeros(1,t_end);
flux_rg  = zeros(1,t_end);
flux_qs  = zeros(1,t_end);
flux_qi  = zeros(1,t_end);
flux_qg  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Tension water
% S2. Free water
% S3. Interflow 
% S4. Baseflow

% RB(aim,P(t)): impervious surface runoff
RB = split_1;

% PI(1-aim,P(t)): precipitation on pervious soil
PI = split_1;

% E(lm,c,S1,Ep(t),delta_t): evaporation from tension water
E = evap_21;

% R(a,b,S1,wmax,PI(1-aim,P(t))): runoff from contributing tension storage
R = saturation_14;

% RS(S2,smax,ex,R(a,b,S1,wmax,PI(1-aim,P(t)))): surface runoff from 
% saturation excess area
RS = saturation_2;

% RI(S2,smax,ex,S2*ki): interflow from contributing area
RI = saturation_2;

% RG(S2,smax,ex,S2*kg): interflow from contributing area
RG = saturation_2;

% QI(ci,S3): interflow
QI = interflow_5;

% QG(cg,S4): baseflow
QG = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0;
                                               1,1,0,0;
                                               0,1,1,0;
                                               0,0,1,1]);                   % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0;
                                                  1,1,0,0;
                                                  0,1,1,0;
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
               (PI(1-aim,P(t)) - ...
                E(lm,c,S1,Ep(t),delta_t) - ...
                R(a,b,S1,wmax,PI(1-aim,P(t))));
            
    tmpf_S2 = @(S1,S2,S3,S4) ...
               (R(a,b,S1,wmax,PI(1-aim,P(t))) - ...
                RS(S2,smax,ex,R(a,b,S1,wmax,PI(1-aim,P(t)))) - ...
                RI(S2,smax,ex,S2*ki) - ...
                RG(S2,smax,ex,S2*kg));
            
    tmpf_S3 = @(S1,S2,S3,S4) ...
               (RI(S2,smax,ex,S2*ki) - ...
                QI(ci,S3));

    tmpf_S4 = @(S1,S2,S3,S4) ... 
               (RG(S2,smax,ex,S2*kg) - ...
                QG(cg,S4));
    
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
                                        [S1old,S2old,S3old,S4old], ...      % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_rb(t) = RB(aim,P(t));
    flux_pi(t) = PI(1-aim,P(t));
    flux_e(t)  = E(lm,c,tmp_sFlux(1),Ep(t),delta_t);
    flux_r(t)  = R(a,b,tmp_sFlux(1),wmax,flux_pi(t));
    flux_rs(t) = RS(tmp_sFlux(2),smax,ex,flux_r(t));
    flux_ri(t) = RI(tmp_sFlux(2),smax,ex,tmp_sFlux(2)*ki);
    flux_rg(t) = RG(tmp_sFlux(2),smax,ex,tmp_sFlux(2)*kg);
    flux_qs(t) = flux_rb(t) + flux_rs(t);
    flux_qi(t) = QI(ci,tmp_sFlux(3));
    flux_qg(t) = QG(cg,tmp_sFlux(4));
    
    % Update the stores
    store_S1(t) = S1old + (flux_pi(t) - flux_e(t)  - flux_r(t)) * delta_t;
    store_S2(t) = S2old + (flux_r(t)  - flux_rs(t) - flux_ri(t) - flux_rg(t)) * delta_t;
    store_S3(t) = S3old + (flux_ri(t) - flux_qi(t)) * delta_t;
    store_S4(t) = S4old + (flux_rg(t) - flux_qg(t)) * delta_t;
  
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_e * delta_t;
    fluxOutput.Q      = (flux_qs + flux_qi + flux_qg) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.rb   = flux_rb * delta_t;
    fluxInternal.pi   = flux_pi * delta_t;
    fluxInternal.r    = flux_r * delta_t;
    fluxInternal.rs   = flux_rs * delta_t;
    fluxInternal.ri   = flux_ri * delta_t;
    fluxInternal.rg   = flux_rg * delta_t;
    fluxInternal.qs   = flux_qs * delta_t;
    fluxInternal.qi   = flux_qi * delta_t;
    fluxInternal.qg   = flux_qg * delta_t;

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






