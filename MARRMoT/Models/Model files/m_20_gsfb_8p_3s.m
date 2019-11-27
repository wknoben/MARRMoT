function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_20_gsfb_8p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [GSFB] 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Nathan, R. J., & McMahon, T. A. (1990). SFB model part l . Validation of 
% fixed model parameters. In Civil Eng. Trans. (pp. 157–161).
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral catchments.
% Water Resources Research, 33(1), 153–166. http://doi.org/doi:10.1029/96WR02840


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
c       = theta(1);     % Recharge time coeffcient [d-1]
ndc     = theta(2);     % Threshold fraction of Smax [-]
smax    = theta(3);     % Maximum soil moisture storage [mm]
emax    = theta(4);     % Maximum evaporation flux [mm/d]
frate   = theta(5);     % Maximum infiltration rate [mm/d]
b       = theta(6);     % Fraction of subsurface flow that is baseflow [-]
dpf     = theta(7);     % Baseflow time coefficient [d-1]
sdrmax  = theta(8);     % Threshold before baseflow can occur [mm]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial near surface storage
S20     = storeInitial(2);       % Initial subsurface storage
S30     = storeInitial(3);       % Initial deepstorage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];             % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);

flux_ea  = zeros(1,t_end);
flux_qs  = zeros(1,t_end);
flux_f   = zeros(1,t_end);
flux_qb  = zeros(1,t_end);
flux_dp  = zeros(1,t_end);
flux_qdr = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Near-surface store
% S2. Subsurface store
% S3. Deep store

% EA(emax,ndc,S1,smax,Ep(t),delta_t): evaporation
EA = evap_20;

% QS(P(t),S1,smax): saturation excess
QS = saturation_1;

% F(frate,ndc*smax,S1,delta_t): infiltration above threshold
F = interflow_11;

% QDR(c,ndc*smax,S3,S1): recharge to upper zone 
QDR = recharge_5;

% QB(b*dpf,sdrmax,S2): baseflow above threshold
QB = baseflow_9;

% DP((1-b)*dpf,S2): deep percolation 
DP = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,1;
                                               1,1,0;
                                               1,1,1]);                       % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,1;
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
    tmpf_S1 = @(S1,S2,S3) ...
            (P(t) + ...
             QDR(c,ndc*smax,S3,S1) - ...
             EA(emax,ndc,S1,smax,Ep(t),delta_t) - ...
             QS(P(t),S1,smax) - ...
             F(frate,ndc*smax,S1,delta_t));
         
    tmpf_S2 = @(S1,S2,S3) ...
            (F(frate,ndc*smax,S1,delta_t) - ...
             QB(b*dpf,sdrmax,S2) - ...
             DP((1-b)*dpf,S2));
         
    tmpf_S3 = @(S1,S2,S3) ...
            (DP((1-b)*dpf,S2) - ...
             QDR(c,ndc*smax,S3,S1));
    
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

    % Calculate the fluxes
    flux_ea(t)  = EA(emax,ndc,tmp_sFlux(1),smax,Ep(t),delta_t);
    flux_qs(t)  = QS(P(t),tmp_sFlux(1),smax);
    flux_f(t)   = F(frate,ndc*smax,tmp_sFlux(1),delta_t);
    flux_qb(t)  = QB(b*dpf,sdrmax,tmp_sFlux(2));
    flux_dp(t)  = DP((1-b)*dpf,tmp_sFlux(2));
    flux_qdr(t) = QDR(c,ndc*smax,tmp_sFlux(3),tmp_sFlux(1));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) + flux_qdr(t) - flux_ea(t) - flux_qs(t) - flux_f(t)) * delta_t;
    store_S2(t) = S2old + (flux_f(t) - flux_qb(t) - flux_dp(t)) * delta_t;
    store_S3(t) = S3old + (flux_dp(t) - flux_qdr(t)) * delta_t;
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = (flux_qs + flux_qb) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.qs   = flux_qs * delta_t;
    fluxInternal.f    = flux_f * delta_t;
    fluxInternal.qb   = flux_qb * delta_t;
    fluxInternal.dp   = flux_dp * delta_t;
    fluxInternal.qdr  = flux_qdr * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...    % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




