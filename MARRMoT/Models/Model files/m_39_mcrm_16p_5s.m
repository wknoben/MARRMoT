function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_39_mcrm_16p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [Midland Catchment Runoff Model] 
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
smax   = theta(1);      % Maximum interception storage [mm]
cmax   = theta(2);      % Maximum fraction of area contributing to rapid runoff [-]
ct     = theta(3);      % Fraction of cmax that is the minimum contributing area c0 [-]
c1     = theta(4);      % Shape parameter for rapid flow distribution [mm-1]
ce     = theta(5);      % Shape parameter for evaporation [mm-1]
dsurp  = theta(6);      % Threshold for direct runoff [mm]
kd     = theta(7);      % Direct runoff time parameter [d-1]
gamd   = theta(8);      % Direct runoff flow non-linearity [-]
qpmax  = theta(9);      % Maximum percolation rate [mm/d]
kg     = theta(10);     % Groundwater time parameter [d-1]
tau    = theta(11);     % Routing delay [d]
sbf    = theta(12);     % Maximum routing store depth [mm]
kcr    = theta(13);     % Channel flow time parameter [d-1]
gamcr  = theta(14);     % Channel flow non-linearity [-]
kor    = theta(15);     % Out-of-bank flow time parameter [d-1]
gamor  = theta(16);     % Out-of-bank flow non-linearity [-]         

% Auxiliary parameters
c0     = ct*cmax;       % Minimum fraction of area contributing to rapid runoff [-]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial interception
S20   = storeInitial(2);       % Initial soil moisture storage
S30   = storeInitial(3);       % Initial groundwater storage
S40   = storeInitial(4);       % Initial channel routing
S50   = storeInitial(5);       % Initial out-of-bank routing

%%DEFINE STORE BOUNDARIES
store_min = [0,-1E6,0,0,0];    % lower bounds of stores
store_upp = [];                % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);

flux_ec  = zeros(1,t_end);
flux_qt  = zeros(1,t_end);
flux_qr  = zeros(1,t_end);
flux_er  = zeros(1,t_end);
flux_qn  = zeros(1,t_end);
flux_qd  = zeros(1,t_end);
flux_qp  = zeros(1,t_end);
flux_qb  = zeros(1,t_end);
flux_uib = zeros(1,t_end);
flux_uob = zeros(1,t_end);
flux_qic = zeros(1,t_end);
flux_qoc = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_full] = uh_7_uniform(1,tau,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_Qt_old  = zeros(1,length(uh_full));   % temporary vector needed to deal with routing

%%INITIALISE FLUX UPDATE VECTOR
tmp_sNew = NaN.*zeros(1,length(storeInitial));

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Soil moisture
% S3. Groundwater
% S4. Channel routing
% S5. Out of bank routing

% EC(S1,Ep(t),delta_t): evaporation from interception
EC = evap_1;

% QT(P(t),S1,smax): overflow from interception
QT = interception_1;

% QR(cmax,c0,c1,S2,QT(P(t),S1,smax)): variable area overland flow
QR = saturation_10;

% QN(QT(P(t),S1,smax),QR(cmax,c0,c1,S2,QT(P(t),S1,smax))): effective 
% infiltration after surface runoff
QN = effective_1;

% ER(ce,S2,Ep(t)- EC(S1,Ep(t),delta_t)): reduced evaporation from soil moisture
ER = evap_17;

% QD(S2,kd,dsurp,gamd,delta_t): non-linear flow above threshold
QD = interflow_9;

% QP(qpmax,dsurp,S2,delta_t): percolation if store is above zero
QP = percolation_6;

% QB(kg,1.5,S3,delta_t): non-linear baseflow
QB = baseflow_7;

% UOB(flux_uib(t),S4,sbf): overflow from channel
UOB = saturation_1;

% QIC(kcr,gamcr,3/4,S4,delta_t): in-channel flow routing
QIC = routing_1;

% QOC(kor,gamor,3/4,S5,delta_t): out-of-bank flow routing
QOC = routing_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
%   Store 1,2 & 3
fsolve_options123    = optimoptions('fsolve','Display','none',...           % [1+n stores] settings of the root finding method
                                     'JacobPattern', [1,0,0;
                                                      1,1,0;
                                                      0,1,1]);              % Specify the Jacobian pattern                                               
lsqnonlin_options123 = optimoptions('lsqnonlin',...                         % lsqnonlin settings for cases when fsolve fails
                                     'Display','none',...
                                     'JacobPattern', [1,0,0;
                                                      1,1,0;
                                                      0,1,1],...
                                     'MaxFunEvals',1000);
%   Store 4 & 5
fsolve_options45     = optimoptions('fsolve','Display','none',...           % [1+n stores] settings of the root finding method
                                     'JacobPattern', [1,0;
                                                      1,1]);                % Specify the Jacobian pattern                                               
lsqnonlin_options45  = optimoptions('lsqnonlin',...                         % lsqnonlin settings for cases when fsolve fails
                                     'Display','none',...
                                     'JacobPattern', [1,0;
                                                      1,1],...
                                     'MaxFunEvals',1000);
                             
%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    
    % NOTE: due to the delay function included in the model, the five
    % stores cannot be solved simultaneously. Instead, stores 1, 2 and 3
    % are solved together to find qr, qd and qb. These fluxes are summed,
    % delayed to form uib and directed towards stores 4 and 5. Then, new
    % storages for stores 4 and 5 can be calculated.

    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    if t == 1; S4old = S40; else; S4old = store_S4(t-1); end                % store 4 at t-1
    if t == 1; S5old = S50; else; S5old = store_S5(t-1); end                % store 5 at t-1

    % ___ Solve stores 1, 2 and 3 _________________________________________
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3)  (P(t) - ...
                            EC(S1,Ep(t),delta_t) - ...
                            QT(P(t),S1,smax));
                        
    tmpf_S2 = @(S1,S2,S3)  (QN(QT(P(t),S1,smax),QR(cmax,c0,c1,S2,QT(P(t),S1,smax))) - ...
                            ER(ce,S2,Ep(t)- EC(S1,Ep(t),delta_t)) - ...
                            QD(S2,kd,dsurp,gamd,delta_t) - ...
                            QP(qpmax,dsurp,S2,delta_t));
    
    tmpf_S3 = @(S1,S2,S3)  (QP(qpmax,dsurp,S2,delta_t) - ...
                            QB(kg,1.5,S3,delta_t));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3);                             % this returns a new anonymous function that we solve in the next step

    % Model solving -------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew(1:3),tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3)),...                  % system of storage equations
                        [S1old,S2old,S3old],...                             % storage values on previous time step
                        fsolve_options123);                                 % solver options
    
    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew(1:3),~,~] = rerunSolver('lsqnonlin', ...                  
                                        lsqnonlin_options123, ...           % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                        eq_sys(1),eq_sys(2),eq_sys(3)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew(1:3), ...                  % recent estimates
                                        [S1old,S2old,S3old], ...            % storages ate previous time step
                                        store_min(1:3), ...                 % lower bounds
                                        store_upp);                         % upper bounds              
    end
    
    % --- Handle the flow routing based on upper store estimates ----------
        % Find the storages needed to update fluxes: update 'tmp_sFlux'
            eval(store_fun);                                                % creates/updates tmp_sFlux 

        % Flow routing
        tmp_Qt_cur    = (QR(cmax,c0,c1,tmp_sFlux(2),QT(P(t),tmp_sFlux(1),smax)) + ...
                         QD(tmp_sFlux(2),kd,dsurp,gamd,delta_t) + ...
                         QB(kg,1.5,tmp_sFlux(3),delta_t)).*uh_full;         % find how the runoff of this time step will be distributed in time
        tmp_Qt_old    = tmp_Qt_old + tmp_Qt_cur;                            % update the 'still-to-flow-out' vector
        flux_uib(t)   = tmp_Qt_old(1);                                      % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_Qt_old    = circshift(tmp_Qt_old,-1);                           % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_Qt_old(end) = 0;                                                % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)

    % --- Solve stores 4 and 5 --------------------------------------------
    % Create temporary store ODE's that need to be solved
    tmpf_S4 = @(S4,S5) (flux_uib(t) - ...
                        UOB(flux_uib(t),S4,sbf) - ...
                        QIC(kcr,gamcr,3/4,S4,delta_t));
    
    tmpf_S5 = @(S4,S5) (UOB(flux_uib(t),S4,sbf) - ...
                        QOC(kor,gamor,3/4,S5,delta_t));

    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S4old,S5old],...
                      delta_t,...
                      tmpf_S4,tmpf_S5);

    % Model solving -------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew(4:5),tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                                    eq_sys(1),eq_sys(2)),...                % system of storage equations
                                    [S4old,S5old],...                       % storage values on previous time step
                                    fsolve_options45);                      % solver options
    
    % --- Check if the solver accuracy ---
    tmp_resnorm = sum(tmp_fval.^2);
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew(4:5),~,~] = rerunSolver('lsqnonlin', ...                  
                                        lsqnonlin_options45, ...            % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                        eq_sys(1),eq_sys(2)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew(4:5), ...                  % recent estimates
                                        [S4old,S5old], ...                  % storages ate previous time step
                                        store_min(4:5), ...                 % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ec(t)  = EC(tmp_sFlux(1),Ep(t),delta_t);
    flux_qt(t)  = QT(P(t),tmp_sFlux(1),smax);
    flux_qr(t)  = QR(cmax,c0,c1,tmp_sFlux(2),flux_qt(t));
    flux_er(t)  = ER(ce,tmp_sFlux(2),Ep(t)- flux_ec(t));
    flux_qn(t)  = QN(flux_qt(t),flux_qr(t));
    flux_qd(t)  = QD(tmp_sFlux(2),kd,dsurp,gamd,delta_t);
    flux_qp(t)  = QP(qpmax,dsurp,tmp_sFlux(2),delta_t);
    flux_qb(t)  = QB(kg,1.5,tmp_sFlux(3),delta_t);
    flux_uob(t) = UOB(flux_uib(t),tmp_sFlux(4),sbf);
    flux_qic(t) = QIC(kcr,gamcr,3/4,tmp_sFlux(4),delta_t);
    flux_qoc(t) = QOC(kor,gamor,3/4,tmp_sFlux(5),delta_t);
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ec(t) - flux_qt(t)) * delta_t;
    store_S2(t) = S2old + (flux_qn(t) - flux_er(t) - flux_qd(t) - flux_qp(t)) * delta_t;
    store_S3(t) = S3old + (flux_qp(t) - flux_qb(t)) * delta_t;
    store_S4(t) = S4old + (flux_uib(t) - flux_uob(t) - flux_qic(t)) * delta_t;
    store_S5(t) = S5old + (flux_uob(t) - flux_qoc(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ec + flux_er) * delta_t;
    fluxOutput.Q      = (flux_qic + flux_qoc) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ec  = flux_ec * delta_t;
    fluxInternal.qt  = flux_qt * delta_t;
    fluxInternal.qr  = flux_qr * delta_t;
    fluxInternal.qn  = flux_qn * delta_t;
    fluxInternal.er  = flux_er * delta_t;
    fluxInternal.qd  = flux_qd * delta_t;
    fluxInternal.qp  = flux_qp * delta_t;
    fluxInternal.qb  = flux_qb * delta_t;
    fluxInternal.uib = flux_uib * delta_t;
    fluxInternal.uob = flux_uob * delta_t;
    fluxInternal.qic = flux_qic * delta_t;
    fluxInternal.qoc = flux_qoc * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     tmp_Qt_old);       % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end



 




