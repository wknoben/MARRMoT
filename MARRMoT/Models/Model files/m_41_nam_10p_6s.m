function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_41_nam_10p_6s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: NAM
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Nielsen, S. A., & Hansen, E. (1973). Numerical simulation of he rainfall-
% runoff process on a daily basis. Nordic Hydrology, (4), 171–190. 
% http://doi.org/https://doi.org/10.2166/nh.1973.0013


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
cs    = theta(1);     % Degree-day factor for snowmelt [mm/oC/d]
cif   = theta(2);     % Runoff coefficient for interflow [d-1]
stot  = theta(3);     % Maximum total soil moisture depth [mm]
cl1   = theta(4);     % Lower zone filling threshold for interflow generation [-]
fl    = theta(5);     % Fraction of total soil depth that makes up lstar
cof   = theta(6);     % Runoff coefficient for overland flow [d-1]
cl2   = theta(7);     % Lower zone filling threshold for overland flow generation [-]
k0    = theta(8);     % Overland flow routing delay [d-1]
k1    = theta(9);     % Interflow routing delay [d-1]
kb    = theta(10);    % Baseflow routing delay [d-1]

% Auxiliary parameters
lstar = fl*stot;      % Maximum lower zone storage [mm]
ustar = (1-fl)*stot;  % Upper zone maximum storage [mm]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial snow pack storage
S20   = storeInitial(2);       % Initial upper soil moisture storage
S30   = storeInitial(3);       % Initial lower soil moisture storage
S40   = storeInitial(4);       % Initial overland flow routing storage
S50   = storeInitial(5);       % Initial interflow routing storage
S60   = storeInitial(6);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0];           % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS 
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);

flux_ps   = zeros(1,t_end);
flux_pr   = zeros(1,t_end);
flux_m    = zeros(1,t_end);
flux_eu   = zeros(1,t_end);
flux_pn   = zeros(1,t_end);
flux_of   = zeros(1,t_end);
flux_inf  = zeros(1,t_end);
flux_if   = zeros(1,t_end);
flux_dl   = zeros(1,t_end);
flux_gw   = zeros(1,t_end);
flux_el   = zeros(1,t_end);
flux_qo   = zeros(1,t_end);
flux_qi   = zeros(1,t_end);
flux_qb   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Snow pack
% S2. Upper soil moisture
% S3. Lower soil moisture
% S4. Overland flow
% S5. Interflow
% S6. Groundwater

% PS(P(t),T(t),0): snowfall
PS = snowfall_1;

% PR(P(t),T(t),0): rainfall
PR = rainfall_1;

% M(cs,0,T(t),S1,delta_t): melt
M = melt_1;

% EU(S2,Ep(t),delta_t): evaporation from upper soil moisture
EU = evap_1;

% PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1),S2,ustar): overflow from upper soil moisture
PN = saturation_1;

% IF(cif,cl1,S2,S3,lstar): interflow dependent on lower soil moisture storage
IF = interflow_6;

% OF(cof,cl2,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1),S2,ustar),S3,lstar): overland flow dependent on lower soil moisture.
OF = interflow_6;

% DL(1-S3/lstar,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1),S2,ustar) - OF(cof,cl2,S2,S3,lstar)):
% fraction flow to lower evaporation zone
DL = split_1;

% GW(S3/lstar,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1),S2,ustar) - OF(cof,cl2,S2,S3,lstar)):
% fraction flow to groundwater
GW = split_1;

% EL(Ep(t),S3,lstar,S2,0.01,delta_t): evaporation from lower zone
EL = evap_15;

% QO(k0,S4): overland flow routing
QO = interflow_5;

% QI(k1,S5): interflow routing
QI = interflow_5;

% QB(kb,S6): baseflow routing
QB = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,1,0,0;
                                               0,1,1,0,1,0;
                                               1,1,1,0,0,1]);               % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,0,0,0,0,0;
                                                    1,1,1,0,0,0;
                                                    1,1,1,0,0,0;
                                                    1,1,1,1,0,0;
                                                    0,1,1,0,1,0;
                                                    1,1,1,0,0,1],...
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
    if t == 1; S6old = S60; else; S6old = store_S6(t-1); end                % store 6 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5,S6) ...
               (PS(P(t),T(t),0) -...
                M(cs,0,T(t),S1,delta_t));

    tmpf_S2 = @(S1,S2,S3,S4,S5,S6) ...
               (PR(P(t),T(t),0) + ...
                M(cs,0,T(t),S1,delta_t) - ...
                EU(S2,Ep(t),delta_t) - ...
                PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar) - ...
                IF(cif,cl1,S2,S3,lstar));
            
	tmpf_S3 = @(S1,S2,S3,S4,S5,S6) ...
               (DL(1-S3/lstar,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar) - OF(cof,cl2,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar),S3,lstar)) - ...
                EL(Ep(t),S3,lstar,S2,0.01,delta_t));
            
    tmpf_S4 = @(S1,S2,S3,S4,S5,S6) ...
               (OF(cof,cl2,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar),S3,lstar) - ...
                QO(k0,S4));
            
    tmpf_S5 = @(S1,S2,S3,S4,S5,S6) ...
               (IF(cif,cl1,S2,S3,lstar) - ...
                QI(k1,S5));
            
    tmpf_S6 = @(S1,S2,S3,S4,S5,S6) ...
               (GW(S3/lstar,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar) - OF(cof,cl2,PN(PR(P(t),T(t),0)+M(cs,0,T(t),S1,delta_t),S2,ustar),S3,lstar)) - ...
                QB(kb,S6));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old,S6old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5,tmpf_S6);     % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),...
                        eq_sys(4),eq_sys(5),eq_sys(6)),...                  % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old,S6old],...           % storage values on previous time step
                        fsolve_options);                                    % solver options

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                        eq_sys(1),eq_sys(2),eq_sys(3),...
                                        eq_sys(4),eq_sys(5),eq_sys(6)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,...
                                         S4old,S5old,S6old], ...            % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
   end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ps(t)   = PS(P(t),T(t),0);
    flux_pr(t)   = PR(P(t),T(t),0);
    flux_m(t)    = M(cs,0,T(t),tmp_sFlux(1),delta_t);
    flux_eu(t)   = EU(tmp_sFlux(2),Ep(t),delta_t);
    flux_pn(t)   = PN(flux_pr(t)+flux_m(t),tmp_sFlux(2),ustar);
    flux_of(t)   = OF(cof,cl2,flux_pn(t),tmp_sFlux(3),lstar);
    flux_inf(t)  = flux_pn(t)-flux_of(t);
    flux_if(t)   = IF(cif,cl1,tmp_sFlux(2),tmp_sFlux(3),lstar);
    flux_dl(t)   = DL(1-tmp_sFlux(3)/lstar,flux_inf(t));
    flux_gw(t)   = GW(tmp_sFlux(3)/lstar,flux_inf(t));
    flux_el(t)   = EL(Ep(t),tmp_sFlux(3),lstar,tmp_sFlux(2),0.01,delta_t);
    flux_qo(t)   = QO(k0,tmp_sFlux(4));
    flux_qi(t)   = QI(k1,tmp_sFlux(5));
    flux_qb(t)   = QB(kb,tmp_sFlux(6));
    
    % Update the stores
    store_S1(t) = S1old + (flux_ps(t) - flux_m(t)) * delta_t;
    store_S2(t) = S2old + (flux_pr(t) + flux_m(t) - flux_eu(t) - flux_if(t) - flux_pn(t)) * delta_t;    
    store_S3(t) = S3old + (flux_dl(t) - flux_el(t)) * delta_t;    
    store_S4(t) = S4old + (flux_of(t) - flux_qo(t)) * delta_t;    
    store_S5(t) = S5old + (flux_if(t) - flux_qi(t)) * delta_t;    
    store_S6(t) = S6old + (flux_gw(t) - flux_qb(t)) * delta_t;    
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eu + flux_el) * delta_t;
    fluxOutput.Q      = (flux_qo + flux_qi + flux_qb) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ps  = flux_ps * delta_t;
    fluxInternal.pr  = flux_pr * delta_t;
    fluxInternal.m   = flux_m * delta_t;
    fluxInternal.eu  = flux_eu * delta_t;
    fluxInternal.pn  = flux_pn * delta_t;
    fluxInternal.of  = flux_of * delta_t;
    fluxInternal.inf = flux_inf * delta_t;
    fluxInternal.if  = flux_if * delta_t;
    fluxInternal.dl  = flux_dl * delta_t;
    fluxInternal.gw  = flux_gw * delta_t;
    fluxInternal.el  = flux_el * delta_t;
    fluxInternal.qo  = flux_qo * delta_t;
    fluxInternal.qi  = flux_qi * delta_t;
    fluxInternal.qb  = flux_qb * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




