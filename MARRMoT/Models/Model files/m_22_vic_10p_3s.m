function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_22_vic_10p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: VIC
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. A., Vrugt, J. A., 
% Gupta, H. V., … Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between 
% hydrological models. Water Resources Research, 44(12), W00B02. 
% http://doi.org/10.1029/2007WR006735
%
% Liang, X., Lettenmaier, D. P., Wood, E. F., & Burges, S. J. (1994). A 
% simple hydrologically based model of land surface water and energy fluxes
% for general circulation models. Journal of Geophysical Research, 99, 
% 14415–14428.


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
ibar    = theta(1);     % Mean interception capacity [mm]
idelta  = theta(2);     % Seasonal interception change as fraction of mean [-]
ishift  = theta(3);     % Maximum interception peak timing [-]
stot    = theta(4);     % Total available storage [mm]
fsm     = theta(5);     % Fraction of stot that constitutes maximum soil mositure storage [-]
b       = theta(6);     % Infiltration excess shape parameter [-]
k1      = theta(7);     % Percolation time parameter [d-1]
c1      = theta(8);     % Percolation non-linearity parameter [-]
k2      = theta(9);     % Baseflow time parameter [d-1]
c2      = theta(10);    % Baseflow non-linearity parameter

% Auxiliary parameter
smmax   = fsm*stot;     % Maximum soil moisture capacity [mm]
gwmax   = (1-fsm)*stot; % Maximum groundwater storage [mm]
tmax    = 365.25;       % Length of one growing cycle [d]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial interception storage
S20     = storeInitial(2);       % Initial soil moisture storage
S30     = storeInitial(3);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];             % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);

flux_ei     = zeros(1,t_end);
flux_peff   = zeros(1,t_end);
flux_iex    = zeros(1,t_end);
flux_qie    = zeros(1,t_end);
flux_inf    = zeros(1,t_end);
flux_et1    = zeros(1,t_end);
flux_qex1   = zeros(1,t_end);
flux_pc     = zeros(1,t_end);
flux_et2    = zeros(1,t_end);
flux_qex2   = zeros(1,t_end);
flux_qb     = zeros(1,t_end);

aux_imax    = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Soil moisture
% S3. Groundwater

% IMAX(ibar,idelta,ishift,t,tmax,delta_t): time-varying maximum interception
IMAX = phenology_2;

% EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t),delta_t): evaporation from interception
EI = evap_7;

% PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)): overflow from interception. 
PEFF = interception_1;

% IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t): excess flow from reduced store capacity.
IEX = excess_1;

% QIE(S2,smmax,b,PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))): 
% infiltration excess flow.
QIE = saturation_2;

% INF(PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)),QIE(S2,smmax,b,PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)))):
% infiltration to soil moisture
INF = effective_1;

% ET1(S2,smmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t))),delta_t):
% evaporation from soil moisture
ET1 = evap_7;

% QEX1(INF(PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)),QIE(S2,smmax,b,PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)))),S2,smmax):
% soil moisture excess
QEX1 = saturation_1;

% PC(k1,c1,S2,smmax,delta_t): percolation from soil moisture
PC = percolation_5;

% ET2(S3,gwmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t))-ET1(S2,smmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t))))),delta_t):
% transpiration from groundwater
ET2 = evap_7;

% QEX2(PC(k1,c1,S2,smmax),S3,gwmax): overflow from groundwater. Threshold
% discontinuity that needs logisitic smoothing.
QEX2 = saturation_1;

% QB(k2,c2,S3,gwmax,delta_t): non-linear baseflow 
QB = baseflow_5;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0;
                                               1,1,0;
                                               1,1,1]);                      % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0;
                                                  1,1,0;
                                                  1,1,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3) ...                                               % store 1 function with current flux values
                (P(t) -...
                 EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t),delta_t) - ...
                 PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t)) - ...
                 IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t));  

    tmpf_S2 = @(S1,S2,S3) ...                                               % store 2 function
                (INF(PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t),QIE(S2,smmax,b,PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t))) - ...
                 ET1(S2,smmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t),delta_t)),delta_t) - ...
                 QEX1(INF(PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t),QIE(S2,smmax,b,PEFF(P(t),S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t))+IEX(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),delta_t))),S2,smmax) - ...
                 PC(k1,c1,S2,smmax,delta_t));                                       

    tmpf_S3 = @(S1,S2,S3) ...                                               % store 3 function
                (PC(k1,c1,S2,smmax,delta_t) - ...
                 ET2(S3,gwmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t),delta_t)-ET1(S2,smmax,max(0,Ep(t)-EI(S1,IMAX(ibar,idelta,ishift,t,tmax,delta_t),Ep(t),delta_t)),delta_t)),delta_t) - ...
                 QEX2(PC(k1,c1,S2,smmax,delta_t),S3,gwmax) - ...
                 QB(k2,c2,S3,gwmax,delta_t));                                   
             
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3);                             % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3)),...                  % system of storage equations
                        [S1old,S2old,S3old],...                             % storage values on previous time step
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
                                        eq_sys(1),eq_sys(2),eq_sys(3)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old], ...            % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % calculate auxiliary variable
    aux_imax(t)    = IMAX(ibar,idelta,ishift,t,tmax,delta_t);

    % Calculate the fluxes
    flux_ei(t)     = EI(tmp_sFlux(1),aux_imax(t),Ep(t),delta_t);
    flux_peff(t)   = PEFF(P(t),tmp_sFlux(1),aux_imax(t));
    flux_iex(t)    = IEX(tmp_sFlux(1),aux_imax(t),delta_t);
    flux_qie(t)    = QIE(tmp_sFlux(2),smmax,b,flux_peff(t)+flux_iex(t));
    flux_inf(t)    = INF(flux_peff(t)+flux_iex(t),flux_qie(t));
    flux_et1(t)    = ET1(tmp_sFlux(2),smmax,max(0,Ep(t)-flux_ei(t)),delta_t);
    flux_qex1(t)   = QEX1(flux_inf(t),tmp_sFlux(2),smmax);
    flux_pc(t)     = PC(k1,c1,tmp_sFlux(2),smmax,delta_t);
    flux_et2(t)    = ET2(tmp_sFlux(3),gwmax,max(0,Ep(t)-flux_ei(t)-flux_et1(t)),delta_t);
    flux_qex2(t)   = QEX2(flux_pc(t),tmp_sFlux(3),gwmax);
    flux_qb(t)     = QB(k2,c2,tmp_sFlux(3),gwmax,delta_t);
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ei(t) - flux_peff(t) - flux_iex(t)) * delta_t;
    store_S2(t) = S2old + (flux_inf(t) - flux_et1(t) - flux_qex1(t) - flux_pc(t)) * delta_t;    
    store_S3(t) = S3old + (flux_pc(t)  - flux_et2(t) - flux_qex2(t) - flux_qb(t)) * delta_t;
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ei + flux_et1 + flux_et2) * delta_t;
    fluxOutput.Q      = (flux_qie + flux_qex1 + flux_qex2 + flux_qb) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ei     = flux_ei * delta_t;
    fluxInternal.et1    = flux_et1 * delta_t;
    fluxInternal.et2    = flux_et2 * delta_t;
    fluxInternal.qie    = flux_qie * delta_t;
    fluxInternal.qex1   = flux_qex1 * delta_t;
    fluxInternal.qex2   = flux_qex2 * delta_t;
    fluxInternal.qb     = flux_qb * delta_t;
    fluxInternal.peff   = flux_peff * delta_t;
    fluxInternal.iex    = flux_iex * delta_t;
    fluxInternal.inf    = flux_inf * delta_t;
    fluxInternal.pc     = flux_pc * delta_t;

    % --- Stores ---
    storeInternal.S1    = store_S1;
    storeInternal.S2    = store_S2;
    storeInternal.S3    = store_S3;
%   storeInternal.imax  = aux_imax;                     % This can be used to return the time-varying maximum interception storage

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




