function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_38_tank2_16p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Tank Model - SMA
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Sugawara, M. (1995). Tank model. In V. P. Singh (Ed.), Computer models of
% watershed hydrology (pp. 165–214). Water Resources Publications, USA.

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
% NOTE1: thresholds for runoff are determine using a single maximum depth
% st for consistency with other models.
% NOTE2: Sugawara specifies that A1 > A2 > B1 > C1 > D1 so that time scales
% increase for deeper stores. This is enforced here by specifying deeper
% time constants (a2,b1,c1,d1) as a fraction of the preceding coefficient.
a0 = theta(1);          % Time parameter for drainage 1>2 [d-1]
b0 = theta(2);          % Time parameter for drainage 2>3 [d-1]
c0 = theta(3);          % Time parameter for drainage 3>4 [d-1]
a1 = theta(4);          % Time parameter for surface runoff 1 [d-1]
fa = theta(5);          % Fraction of a1 that is a2 [-]
fb = theta(6);          % Fraction of a2 that is b1 [-]   
fc = theta(7);          % Fraction of b1 that is c1 [-]
fd = theta(8);          % Fraction of c1 that is d1 [-]
st = theta(9);          % Maximum soil depth (sum of runoff thresholds) [mm]
f2 = theta(10);         % Fraction of st-sm1 that is added to sm1 to find threshold t2 [-] (ensures t2 > sm1)
f1 = theta(11);         % Fraction of st-t2 that is added to t2 to find threshold 1 [-] (ensures t1 > t2)
f3 = theta(12);         % Fraction of st-t1-sm2 that consitutes threshold 3 [-]
k1 = theta(13);         % Base rate of capillary rise [mm/d]
k2 = theta(14);         % Base rate of soil moisture exchange [mm/d]
z1 = theta(15);         % Fraction st that is sm1 [-]
z2 = theta(16);         % Fraction of st-t1 that is sm2 [-]

% Auxiliary parameters
sm1 = z1*st;            % Size of primary soil moisture store, threshold before F12 starts [mm]
t2  = sm1+f2*(st-sm1);  % Threshold before surface runoff Y2 starts [mm]
t1  = t2+f1*(st-t2);    % Thresold before surface runoff Y1 starts [mm]
sm2 = z2*(st-t1);       % Size of secondary soil moisture store S5 [mm]
t3  = f3*(st-t1-sm2);   % Threshold before intermediate runoff starts [mm]
t4  = st-t1-sm2-t3;     % Threshold before sub-base runoff starts [mm]
a2  = fa*a1;            % Time parameter for surface runoff 2 [d-1]
b1  = fb*a2;            % Time parameter for intermediate runoff 1 [d-1]
c1  = fc*b1;            % Time parameter for sub-base runoff 1 [d-1]
d1  = fd*c1;            % Time parameter for base runoff 1 [d-1]

%%INITIALISE MODEL STORES
S10 = storeInitial(1);  % Initial upper layer storage
S20 = storeInitial(2);  % Initial second layer storage
S30 = storeInitial(3);  % Initial third layer storage
S40 = storeInitial(4);  % Initial fourth layer storage
S50 = storeInitial(5);  % Initial secondary soil moisture layer storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0];             % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);

flux_t1  = zeros(1,t_end);
flux_t2  = zeros(1,t_end);
flux_y1  = zeros(1,t_end);
flux_y2  = zeros(1,t_end);
flux_y3  = zeros(1,t_end);
flux_y4  = zeros(1,t_end);
flux_y5  = zeros(1,t_end);
flux_e1  = zeros(1,t_end);
flux_e2  = zeros(1,t_end);
flux_e3  = zeros(1,t_end);
flux_e4  = zeros(1,t_end);
flux_f12 = zeros(1,t_end);
flux_f23 = zeros(1,t_end);
flux_f34 = zeros(1,t_end);


%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Upper layer storage
% S2. Second layer storage
% S3. Third layer storage
% S4. Fourth layer storage
% S5. Secondary soil moisture storage

% T1(k1,sm1,S1,S2,delta_t): capillary rise from S2
T1 = capillary_3;

% T2(k2,S1,sm1,S5,sm2): exchange between S1 and S5
T2 = exchange_2;

% E1(S1,Ep(t),delta_t): evaporation from store 1
E1 = evap_1;

% F12(S1,a0,sm1): drainage to layer 2
F12 = interflow_8;

% Y2(S1,a2,t2): threshold excess flow
Y2 = interflow_8;

% Y1(S1,a1,t1): threshold excess flow
Y1 = interflow_8;

% E2(S2,max(0,Ep(t)- E1(S1,Ep(t)),delta_t): evaporation from store 2
E2 = evap_1;

% F23(b0,S2): drainage to layer 3
F23 = recharge_3;

% Y3(S2,b1,t3): threshold excess flow
Y3 = interflow_8;

% E3(S3,max(0,Ep(t)- E1(S1,Ep(t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t)))),delta_t): 
% evaporation from store 3
E3 = evap_1;

% F34(c0,S3): drainage to layer 4.
F34 = recharge_3;

% Y4(S3,c1,t4): threshold excess flow
Y4 = interflow_8;

% E4(S4,max(0,Ep(t)- E1(S1,Ep(t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t)))-E3(S3,max(0,Ep(t)- E1(S1,Ep(t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t)))))),delta_t): 
% evaporation from store 4
E4 = evap_1;

% Y5(d1,S4): baseflow from lower layer
Y5 = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,1;
                                               1,1,0,0,0;
                                               1,1,1,0,0;
                                               1,1,1,1,0;
                                               1,0,0,0,1]);                 % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,0,0,0,1;
                                                    1,1,0,0,0;
                                                    1,1,1,0,0;
                                                    1,1,1,1,0;
                                                    1,0,0,0,1],...
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
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5) ...
                (P(t) + ...
                 T1(k1,sm1,S1,S2,delta_t) - ...
                 T2(k2,S1,sm1,S5,sm2) - ...
                 E1(S1,Ep(t),delta_t) - ...
                 F12(S1,a0,sm1) - ...
                 Y2(S1,a2,t2) - ...
                 Y1(S1,a1,t1));                                             % Store 1 function with current flux values
    
    tmpf_S2 = @(S1,S2,S3,S4,S5) ...
                (F12(S1,a0,sm1) - ...
                 T1(k1,sm1,S1,S2,delta_t) - ...
                 E2(S2,max(0,Ep(t)- E1(S1,Ep(t),delta_t)),delta_t) - ...
                 F23(b0,S2) - ...
                 Y3(S2,b1,t3));                                             % Store 2 function
    
    tmpf_S3 = @(S1,S2,S3,S4,S5) ...
                (F23(b0,S2) - ...
                 E3(S3,max(0,Ep(t)- E1(S1,Ep(t),delta_t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t),delta_t)),delta_t)),delta_t) - ...
                 F34(c0,S3) - ...
                 Y4(S3,c1,t4));                                             % Store 3 function
             
    tmpf_S4 = @(S1,S2,S3,S4,S5) ...
                (F34(c0,S3) - ...
                 E4(S4,max(0,Ep(t)-E1(S1,Ep(t),delta_t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t),delta_t)),delta_t)-E3(S3,max(0,Ep(t)- E1(S1,Ep(t),delta_t)-E2(S2,max(0,Ep(t)- E1(S1,Ep(t),delta_t)),delta_t)),delta_t)),delta_t) - ...
                 Y5(d1,S4));                                                % Store 4 function
    
    tmpf_S5 = @(S1,S2,S3,S4,S5) ...
                (T2(k2,S1,sm1,S5,sm2));                                     % Store 5 function
             
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5);             % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                            eq_sys(5)),...                                  % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old],...                 % storage values on previous time step
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
                                        eq_sys(1),eq_sys(2),eq_sys(3),...
                                            eq_sys(4),eq_sys(5)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old,S5old], ...% storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_t1(t)  = T1(k1,sm1,tmp_sFlux(1),tmp_sFlux(2),delta_t);
    flux_t2(t)  = T2(k2,tmp_sFlux(1),sm1,tmp_sFlux(5),sm2);
    flux_y1(t)  = Y1(tmp_sFlux(1),a1,t1);
    flux_y2(t)  = Y2(tmp_sFlux(1),a2,t2);
    flux_y3(t)  = Y3(tmp_sFlux(2),b1,t3);
    flux_y4(t)  = Y4(tmp_sFlux(3),c1,t4);
    flux_y5(t)  = Y5(d1,tmp_sFlux(4));
    flux_e1(t)  = E1(tmp_sFlux(1),Ep(t),delta_t);
    flux_e2(t)  = E2(tmp_sFlux(2),max(0,Ep(t)-flux_e1(t)),delta_t);
    flux_e3(t)  = E3(tmp_sFlux(3),max(0,Ep(t)-flux_e1(t)-flux_e2(t)),delta_t);
    flux_e4(t)  = E4(tmp_sFlux(4),max(0,Ep(t)-flux_e1(t)-flux_e2(t)-flux_e3(t)),delta_t);
    flux_f12(t) = F12(tmp_sFlux(1),a0,sm1);
    flux_f23(t) = F23(b0,tmp_sFlux(2));
    flux_f34(t) = F34(c0,tmp_sFlux(3));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) + flux_t1(t) - flux_t2(t) - flux_e1(t) - flux_f12(t) - flux_y1(t) - flux_y2(t)) * delta_t;
    store_S2(t) = S2old + (flux_f12(t) - flux_t1(t) - flux_e2(t) - flux_f23(t) - flux_y3(t)) * delta_t;    
    store_S3(t) = S3old + (flux_f23(t) - flux_e3(t) - flux_f34(t) - flux_y4(t)) * delta_t;
    store_S4(t) = S4old + (flux_f34(t) - flux_e4(t) - flux_y5(t)) * delta_t;
    store_S5(t) = S5old + (flux_t2(t)) * delta_t;
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_e1+flux_e2+flux_e3+flux_e4) * delta_t;
    fluxOutput.Q      = (flux_y1+flux_y2+flux_y3+flux_y4+flux_y5) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.t1   = flux_t1 * delta_t;
    fluxInternal.t2   = flux_t2 * delta_t;
    fluxInternal.y1   = flux_y1 * delta_t;
    fluxInternal.y2   = flux_y2 * delta_t;
    fluxInternal.y3   = flux_y3 * delta_t;
    fluxInternal.y4   = flux_y4 * delta_t;
    fluxInternal.y5   = flux_y5 * delta_t;
    fluxInternal.e1   = flux_e1 * delta_t;
    fluxInternal.e2   = flux_e2 * delta_t;
    fluxInternal.e3   = flux_e3 * delta_t;
    fluxInternal.e4   = flux_e4 * delta_t;
    fluxInternal.f12  = flux_f12 * delta_t;
    fluxInternal.f23  = flux_f23 * delta_t;
    fluxInternal.f34  = flux_f34 * delta_t;

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
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




