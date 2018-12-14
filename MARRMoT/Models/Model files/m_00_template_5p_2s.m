function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
    m_00_template_5p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [xxx] 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
%   [reference]


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
kc      = theta(2);     % Maximum capillary rise [mm/d]
kp      = theta(3);     % Maximum percolation [mm/d]
ks      = theta(4);     % Runoff coefficient [d-1]
delay   = theta(5);     % Routing delay [d]
% ...

%%INITIALISE MODEL STORES
S10         = storeInitial(1);       % Initial soil moisture storage
S20         = storeInitial(2);       % Initial groundwater storage
% ...

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                   % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
% ...

flux_cap  = zeros(1,t_end);
flux_ea   = zeros(1,t_end);
flux_qo   = zeros(1,t_end);
flux_perc = zeros(1,t_end);
flux_qs   = zeros(1,t_end);
flux_qt   = zeros(1,t_end);
% ...

%%PREPARE UNIT HYDROGRAPHS
% [Optional]
[~,uh_full] = uh_4_full(1,delay,delta_t);
% ...

%%INITIALISE ROUTING VECTORS
tmp_Qt_old  = zeros(1,length(uh_full));  % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture
% S2. Groundwater

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

% Ea: evaporation from soil moisture. Angle discontinuity
EA = @(S1,Ep,delta_t) min(S1/delta_t,Ep);

% Qo: overflow from soil moisture. This formula uses a threshold - this 
% gives a threshold discontinuity which we deal with using a logistic smoother
QO = @(P,S1,S1max) P.*(1-smoothThreshold_storage_logistic(S1,S1max,0.001));

% cap: capillary rise frm groundwater to soil moisture. This can use min 
% function, this leads to an angle discontinuity
CAP = @(kc,S1,S1max,S2,delta_t) min(max(kc*(S1max-S1)/S1max,0),S2/delta_t);

% perc: percolation from soil moisture to groundwater. Angle discontinuity at S1 = 0
PERC = @(kp,S1,S1max,delta_t) min(kp.*S1/S1max,S1/delta_t);

% Qs: flow from groundwater. Angle discontinuity at S2 = 0
QS = @(ks,S2) ks*S2;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1;
                                               1,1]);                       % Specify the Jacobian pattern                                               
% fzero_options = optimset('Display','off');                                % [1 store] settings of the root finding method
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,1;
                                                  1,1],...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = ...
        @(S1,S2) ...                          % Change in S1 depends on ...
         (P(t) + ...                          % Precipitation             +
          CAP(kc,S1,S1max,S2,delta_t) - ...   % Capillary rise to S1      -
          EA(S1,Ep(t),delta_t) - ...          % Evaporation from S1       -
          QO(P(t),S1,S1max) - ...             % Surface runoff from S1    -
          PERC(kp,S1,S1max,delta_t));         % Percolation from S1

    tmpf_S2 = ...
        @(S1,S2) ...                          % Change in S2 depends on ...
         (PERC(kp,S1,S1max,delta_t) - ...     % Percolation to S2         -
          QS(ks,S2) - ...                     % Slow flow from S2         -
          CAP(kc,S1,S1max,S2,delta_t));       % Capillary rise from S2
    
    % Call the numerical scheme function to create the ODE approximations.
    % This returns a new anonymous function that we solve in the next step.
    solve_fun = feval(scheme,...              % numerical approximation method
                      [S1old,S2old],...       % Store values at t-1
                      delta_t,...             % time step size
                      tmpf_S1,tmpf_S2);       % anonymous functions of ODEs    

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...                    % call solve_fun as anonymous function
                        eq_sys(1),eq_sys(2)),...                            % system of storage equations
                        [S1old,S2old],...                                   % storage values on previous time step
                        fsolve_options);                                    % solver options
    
%     [tmp_sNew, tmp_fval] = fzero(solve_fun,...
%                                  S1old,...
%                                  fzero_options);                          % 1 store solver

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       
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
    % This line creates/updates a variable called 'tmp_sFlux' which is used
    % to update the model fluxes for the current time step. Which variables
    % get assigned to 'tmp_sFlux' is a feature of the chosen numerical time
    % stepping scheme (see line 123-124).
    eval(store_fun);   

    % Calculate the fluxes
    flux_cap(t)  = CAP(kc,tmp_sFlux(1),S1max,tmp_sFlux(2),delta_t);
    flux_ea(t)   = EA(tmp_sFlux(1),Ep(t),delta_t);
    flux_qo(t)   = QO(P(t),tmp_sFlux(1),S1max);
    flux_perc(t) = PERC(kp,tmp_sFlux(1),S1max,delta_t);
    flux_qs(t)   = QS(ks,tmp_sFlux(2));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) + flux_cap(t) - flux_ea(t) - ...
                            flux_qo(t) - flux_perc(t)) * delta_t;
    store_S2(t) = S2old + (flux_perc(t) - flux_qs(t) - ...
                            flux_cap(t)) * delta_t;    

% Routing -----------------------------------------------------------------    
    % Total runoff Qt = Qo + Qs. Apply a triangular routing scheme with
    % time base 'delay' (parameter 5)
    tmp_Qt_cur    = (flux_qo(t) + flux_qs(t)).*uh_full;                     % find how the runoff of this time step will be distributed in time
    tmp_Qt_old    = tmp_Qt_old + tmp_Qt_cur;                                % update the 'still-to-flow-out' vector
    flux_qt(t)    = tmp_Qt_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_Qt_old = circshift(tmp_Qt_old,-1);                                  % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_Qt_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the function and should NOT be renamed
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.cap  = flux_cap  * delta_t;
    fluxInternal.perc = flux_perc * delta_t;
    fluxInternal.Qo   = flux_qo   * delta_t;
    fluxInternal.Qs   = flux_qs   * delta_t;

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
                                                        
    % Code for manual water balance calculation
%     waterBalance = sum(P) - ...                         % Precipitation
%                    sum(fluxOutput.Ea) - ...             % Evaporation
%                    sum(fluxOutput.Q) - ...              % Outflow
%                    (store_S1(end)-S10) - ...            % Storage change
%                    (store_S2(end)-S20) - ...
%                    sum(tmp_Qt_old);                     % Still in routing vector
    
%     disp(['Total P  = ',num2str(sum(P)),'mm.'])
%     disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
%     disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
%     disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
%     disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
%     disp(['On route = ',num2str(sum(tmp_Qt_old)),'mm.'])
%     disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(routing)) - delta S = ',num2str(waterBalance)]);  
end




 




