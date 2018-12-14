function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
    m_nn_example_7p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [MARRMoT User Manual example model] 
%   
% Model reference
% MARRMoT User Manual, 2018. 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

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
crate = theta(1);     % Maximum capillary rise rate [mm/d]
uzmax = theta(2);     % Maximum upper zone storage [mm]
prate = theta(3);     % Maximum percolation rate [mm/d]
klz   = theta(4);     % Lower zone runoff coefficient [d-1]
alpha = theta(5);     % Fraction of lower zone runoff to groundwater [-]
kg    = theta(6);     % Groundwater runoff coefficient [d-1]
d     = theta(7);     % Routing delay [d]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial upper zone storage
S20   = storeInitial(2);       % Initial lower zone storage
S30   = storeInitial(3);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];           % lower bounds of stores
store_upp = [];                % optional higher bounds

%%INITIALISE STORAGE VECTORS 
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);

flux_qse = zeros(1,t_end);
flux_e   = zeros(1,t_end);
flux_qp  = zeros(1,t_end);
flux_qc  = zeros(1,t_end);
flux_qlz = zeros(1,t_end);
flux_qf  = zeros(1,t_end);
flux_qg  = zeros(1,t_end);
flux_qs  = zeros(1,t_end);
flux_qt  = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_full] = uh_4_full(1,d,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_Qt_old  = zeros(1,length(uh_full));      % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Upper zone
% S2. Lower zone
% S3. Groundwater

% Model smoothing
% With Matlab's fsolve, smoothing is only needed when the function is 
% undefined (i.e. has thresholds). Angle discontinuities (such as from 
% min(0,x)) can be dealt with by the solver. Thus, threshold
% discontinuities are smoothed with a logisitic function (e.g. Kavetski and
% Kuczera, 2007) with default smoothing parameters (Clark et al, 2008).
%
% Kavetski and Kuczera, 2007. Model smoothing strategies to remove
% microscale discontinuities and spurious secondary optima in objective
% functions in hydrological calibration. Water Resources Research, 43,
% W03411, doi:10.1029/2006WR005195.
%
% Clark, Slater, Rupp, Woods, Vrugt, Gupta, Wagener and Hay, 2008. 
% Framework for Understanding Structural Errors (FUSE): A modular framework
% to diagnose differences between hydrological models. Water Resources 
% Research, 44, doi:10.1029/2007/WR006735.

% E(S1,uzmax,Ep(t),delta_t): evaporation from upper zone (S1). 
E = evap_7;

% QSE(P(t),S1,uzmax): saturation excess from upper zone (S1). 
% Has a threshold discontinuity and needs logistic smoothing
QSE = saturation_1;

% QP(prate,S1,delta_t): percolation from upper zone (S1) to lower zone (S2)
QP = percolation_1;

% QC(crate,S1,uzmax,S2,delta_t): capillary rise from lower (S2) to upper 
% zone (S1)
QC = capillary_1;

% QLZ(klz,S2): outflow from lower zone (S2)
QLZ = baseflow_1;

% QF(1-alpha,QLZ(klz,S2)): fraction (1-alpha) of lower zone outflow (QLZ) 
% that is fast flow
QF = split_1;

% QG(alpha,QLZ(klz,S2)): fraction (alpha) of lower zone outflow (QLZ) that 
% goes to groundwater (S3)
QG = split_1;

% QS(kg,S3): outflow from groundwater (S3)
QS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1,0;
                                               1,1,0;
                                               0,1,1]);                     % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,1,0;
                                                  1,1,0;
                                                  0,1,1],...
                                 'MaxFunEvals',1000);

% Prepare the options for the solver (saves time later)
[fsolve_options,optionFeedback] = ...
    prepareOptionsForSolver(fsolve_options, 'fsolve');     

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = ...
        @(S1,S2,S3) ...                       % Change in S1 depends on ...
         (P(t) + ...                          % Precipitation to S1      +
          QC(crate,S1,uzmax,S2,delta_t) - ... % Capillary rise to S1     -
          E(S1,uzmax,Ep(t),delta_t) - ...     % Evaporation from S1      -
          QSE(P(t),S1,uzmax) - ...            % Surface runoff from S1   -
          QP(prate,S1,delta_t));              % Percolation from S1        

    tmpf_S2 = ...
        @(S1,S2,S3) ...                       % Change in S2 depends on ...
         (QP(prate,S1,delta_t) - ...          % Percolation to S2        -     
          QC(crate,S1,uzmax,S2,delta_t) - ... % Capillary rise from S2   -
          QLZ(klz,S2));                       % Lower zone outflow from S2       

    tmpf_S3 = ...
        @(S1,S2,S3) ...                       % Change in S2 depends on ...
         (QG(alpha,QLZ(klz,S2)) - ...         % Recharge to S3           -
          QS(kg,S3));                         % Slow flow from S3
    
    % Call the numerical scheme function to create the ODE approximations.
    % This returns a new anonymous function that we solve in the next step.
    solve_fun = feval(scheme,...                % numerical approximation method
                      [S1old,S2old,S3old],...   % Store values at t-1
                      delta_t,...               % time step size
                      tmpf_S1,tmpf_S2,tmpf_S3); % anonymous functions of ODEs

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve_noMSG(@(eq_sys) solve_fun(...              % call solve_fun as anonymous function
                        eq_sys(1),eq_sys(2),eq_sys(3)),...                  % system of storage equations
                        [S1old,S2old,S3old],...                             % storage values on previous time step
                        fsolve_options,optionFeedback);                     % solver options
    
    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                          eq_sys(1),eq_sys(2),...
                                          eq_sys(3)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old], ...            % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % This line creates/updates a variable called 'tmp_sFlux' which is used
    % to update the model fluxes for the current time step. Which variables
    % get assigned to 'tmp_sFlux' is a feature of the chosen numerical time
    % stepping scheme (see line 133-134).
    eval(store_fun);                                                        

    % Calculate the fluxes
    flux_qse(t) = QSE(P(t),tmp_sFlux(1),uzmax);
    flux_e(t)   = E(tmp_sFlux(1),uzmax,Ep(t),delta_t);
    flux_qp(t)  = QP(prate,tmp_sFlux(1),delta_t);
    flux_qc(t)  = QC(crate,tmp_sFlux(1),uzmax,tmp_sFlux(2),delta_t);
    flux_qlz(t) = QLZ(klz,tmp_sFlux(2));
    flux_qf(t)  = QF(1-alpha,flux_qlz(t));
    flux_qg(t)  = QG(alpha,flux_qlz(t));
    flux_qs(t)  = QS(kg,tmp_sFlux(3));
        
    % Update the stores
    store_S1(t) = S1old + (P(t)       + flux_qc(t) - flux_e(t) - ...
                            flux_qse(t) - flux_qp(t)) * delta_t;
    store_S2(t) = S2old + (flux_qp(t) - flux_qc(t) - ...
                            flux_qlz(t)) * delta_t;  
    store_S3(t) = S3old + (flux_qg(t) - flux_qs(t)) * delta_t;

% Routing -----------------------------------------------------------------    
    % Total runoff Q = Qse + Qf + Qs. Apply a pre-determined (line 82)
    % triangular Unit Hydrograph routing scheme to find lagged flow Qt.
    tmp_Qt_cur      = (flux_qse(t) + flux_qf(t) + flux_qs(t)).*uh_full;     % find how the runoff of this time step will be distributed in time
    tmp_Qt_old      = tmp_Qt_old + tmp_Qt_cur;                              % update the 'still-to-flow-out' vector
    flux_qt(t)      = tmp_Qt_old(1);                                        % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_Qt_old      = circshift(tmp_Qt_old,-1);                             % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_Qt_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the 
    % funcion and should NOT be renamed
    fluxOutput.Ea     = flux_e  * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qp   = flux_qp  * delta_t;
    fluxInternal.qc   = flux_qc  * delta_t;
    fluxInternal.qlz  = flux_qlz * delta_t;
    fluxInternal.qf   = flux_qf  * delta_t;
    fluxInternal.qg   = flux_qg  * delta_t;
    fluxInternal.qs   = flux_qs  * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;

% Check water balance
if nargout == 4
    waterBalance = ...
     checkWaterBalance(...
      P,...              % Incoming precipitation
      fluxOutput,...     % Fluxes Q and Ea leaving the model
      storeInternal,...  % Time series of storages ...
      storeInitial,...   % And initial store values to calculate delta S
      tmp_Qt_old);       % Whether the model uses a routing scheme that
                         % still contains water. Use '0' for no routing
end


 




