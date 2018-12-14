function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_40_smar_8p_6s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: SMAR
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% O’Connell, P. E., Nash, J. E., & Farrell, J. P. (1970). River flow 
% forecasting through conceptual models part II - the Brosna catchment at 
% Ferbane. Journal of Hydrology, 10, 317–329.
%
% Tan, B. Q., & O’Connor, K. M. (1996). Application of an empirical 
% infiltration equation in the SMAR conceptual model. Journal of Hydrology,
% 185(1-4), 275–295. http://doi.org/10.1016/0022-1694(95)02993-1

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
h    = theta(1);     % Maximum fraction of direct runoff [-] 
y    = theta(2);     % Infiltration rate [mm/d] 
smax = theta(3);     % Maximum soil moisture storage [mm]
c    = theta(4);     % Evaporation reduction coefficient [-]
g    = theta(5);     % Groundwater recharge coefficient [-]
kg   = theta(6);     % Groundwater time parameter [d-1]
n    = theta(7);     % Number of Nash cascade reservoirs [-]
nk   = theta(8);     % Routing delay [d]. n and k are optimized together 
                     % due to interactions between them
k    = nk/n;         % time parameter in the gamma function                     

%%INITIALISE MODEL STORES
S10  = storeInitial(1);       % Initial soil moisture storage 1
S20  = storeInitial(2);       % Initial soil moisture storage 2
S30  = storeInitial(3);       % Initial soil moisture storage 3
S40  = storeInitial(4);       % Initial soil moisture storage 4
S50  = storeInitial(5);       % Initial soil moisture storage 5
S60  = storeInitial(6);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0];           % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1    = zeros(1,t_end);
store_S2    = zeros(1,t_end);
store_S3    = zeros(1,t_end);
store_S4    = zeros(1,t_end);
store_S5    = zeros(1,t_end);
store_S6 	= zeros(1,t_end);

flux_pstar  = zeros(1,t_end);
flux_estar  = zeros(1,t_end);
flux_evap   = zeros(1,t_end);
flux_r1     = zeros(1,t_end);
flux_r2     = zeros(1,t_end);
flux_i      = zeros(1,t_end);
flux_e1     = zeros(1,t_end);
flux_e2     = zeros(1,t_end);
flux_e3     = zeros(1,t_end);
flux_e4     = zeros(1,t_end);
flux_e5     = zeros(1,t_end);
flux_q1     = zeros(1,t_end);
flux_q2     = zeros(1,t_end);
flux_q3     = zeros(1,t_end);
flux_q4     = zeros(1,t_end);
flux_r3     = zeros(1,t_end);
flux_rg     = zeros(1,t_end);
flux_r3star = zeros(1,t_end);
flux_qr     = zeros(1,t_end);
flux_qg     = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_gamma,frac_routing_beyond_time_series] = uh_6_gamma(1,n,k,t_end,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_Qr_old  = zeros(1,length(uh_gamma));            % temporary vector needed to deal with routing
tmp_Qr_beyond_tEnd = zeros(1,length(uh_gamma));     % temporary vector needed track routing beyond simulation time end

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture 1
% S2. Soil moisture 2
% S3. Soil moisture 3
% S4. Soil moisture 4
% S5. Soil moisture 5
% S6. Groundwater

% PSTAR(P(t),Ep(t)): effective rainfall
PSTAR = effective_1;

% ESTAR(Ep(t),P(t)): effective PET
ESTAR = effective_1;

% R1(h,(S1+S2+S3+S4+S5),smax,PSTAR(P(t),ESTAR(Ep(t),P(t)))): excess runoff. 
R1 = saturation_6; 

% I(PSTAR(P(t),ESTAR(Ep(t),P(t)))-R1(h,(S1+S2+S3+S4+S5),smax,PSTAR(P(t),ESTAR(Ep(t),P(t)))) ,y): 
% infiltration into soil moisture:
I = infiltration_4; 

% R2(flux_pstar(t)-flux_r1(t),flux_i(t)): infiltration excess. Used only
% in the flux and routing section
R2 = effective_1;

% E1(c,0,ESTAR(Ep(t),P(t)),S1,delta_t): evaporation from S1
E1 = evap_13;

% Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5):
% overflow to S2
Q1 = saturation_1;

% E2(c,1,ESTAR(Ep(t),P(t)),S2,S1,0.1,delta_t): evaporation from S2
E2 = evap_14;

% Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5);
% overflow from soil moisture 2
Q2 = saturation_1;

% E3(c,2,ESTAR(Ep(t),P(t)),S3,S2,0.1,delta_t): evaporation from S3
E3 = evap_14;

% Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5):
% overflow from soil moisture 3
Q3 = saturation_1;

% E4(c,3,ESTAR(Ep(t),P(t)),S4,S3,0.1,delta_t): evaporation from S4
E4 = evap_14;

% Q4(Q3(Q2(Q1(I(PSTAR(P(t),ESTAR(Ep(t),P(t)))-R1(h,(S1+S2+S3+S4+S5),smax,PSTAR(P(t)ESTAR(Ep(t),P(t)))),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5):
% overflow from soil moisture 4
Q4 = saturation_1;

% E5(c,4,ESTAR(Ep(t),P(t)),S5,S4,0.1,delta_t): evaporation from S5
E5 = evap_14;

% R3(Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5),S5,smax/5):
% overflow from soil moisture 5
R3 = saturation_1;

% RG(g,R3(Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5),S5,smax/5)):
% flow into groundwater reservoir
RG = split_1;

% R3STAR(1-g,R3(Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5),S5,smax/5)):
% flow towards routing scheme
R3STAR = split_1;

% QG(kg,S6)
QG = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0;
                                               1,1,0,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,1,0,0;
                                               1,1,1,1,1,0;
                                               1,1,1,1,1,1]);               % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0,0,0;
                                               1,1,0,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,1,0,0;
                                               1,1,1,1,1,0;
                                               1,1,1,1,1,1],...
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
        
    % Calculate the preliminary fluxes
    flux_pstar(t)  = PSTAR(P(t),Ep(t));
    flux_estar(t)  = ESTAR(Ep(t),P(t));
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5,S6) ...
               (I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y) - ...
                E1(c,0,flux_estar(t),S1,delta_t) - ...
                Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5));                                 % store 1 function with current flux values
    tmpf_S2 = @(S1,S2,S3,S4,S5,S6) ...
               (Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5) - ...
                E2(c,1,flux_estar(t),S2,S1,0.1,delta_t) - ...
                Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5));                               % store 2 function
    tmpf_S3 = @(S1,S2,S3,S4,S5,S6) ...
               (Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5) - ...
                E3(c,2,flux_estar(t),S3,S2,0.1,delta_t) - ...
                Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5));
    tmpf_S4 = @(S1,S2,S3,S4,S5,S6) ...
                (Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5) - ...
                 E4(c,3,flux_estar(t),S4,S3,0.1,delta_t) - ...
                 Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5));
    tmpf_S5 = @(S1,S2,S3,S4,S5,S6) ...
                (Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5) - ...
                 E5(c,4,flux_estar(t),S5,S4,0.1,delta_t) - ...
                 R3(Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5),S5,smax/5));
    tmpf_S6 = @(S1,S2,S3,S4,S5,S6) ...
                (RG(g,R3(Q4(Q3(Q2(Q1(I(flux_pstar(t)-R1(h,(S1+S2+S3+S4+S5),smax,flux_pstar(t)),y),S1,smax/5),S2,smax/5),S3,smax/5),S4,smax/5),S5,smax/5)) - ...
                 QG(kg,S6));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old,S6old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5,tmpf_S6);     % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) ...
                            solve_fun(eq_sys(1),eq_sys(2),...
                                      eq_sys(3),eq_sys(4),...
                                      eq_sys(5),eq_sys(6)),...              % system of storage equations
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
                                                    eq_sys(1),eq_sys(2),...
                                                    eq_sys(3),eq_sys(4),...
                                                    eq_sys(5),eq_sys(6)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,...
                                         S4old,S5old,S6old], ...            % storages at previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_evap(t)   = min(Ep(t),P(t));
    flux_r1(t)     = R1(h,(tmp_sFlux(1)+tmp_sFlux(2)+tmp_sFlux(3)+tmp_sFlux(4)+tmp_sFlux(5)),smax,flux_pstar(t));
    flux_i(t)      = I(flux_pstar(t)-flux_r1(t),y);
    flux_r2(t)     = R2(flux_pstar(t)-flux_r1(t),flux_i(t));
    flux_e1(t)     = E1(c,0,flux_estar(t),tmp_sFlux(1),delta_t);
    flux_e2(t)     = E2(c,1,flux_estar(t),tmp_sFlux(2),tmp_sFlux(1),0.1,delta_t);
    flux_e3(t)     = E3(c,2,flux_estar(t),tmp_sFlux(3),tmp_sFlux(2),0.1,delta_t);
    flux_e4(t)     = E4(c,3,flux_estar(t),tmp_sFlux(4),tmp_sFlux(3),0.1,delta_t);
    flux_e5(t)     = E5(c,4,flux_estar(t),tmp_sFlux(5),tmp_sFlux(4),0.1,delta_t);
    flux_q1(t)     = Q1(flux_i(t), tmp_sFlux(1),smax/5);
    flux_q2(t)     = Q2(flux_q1(t),tmp_sFlux(2),smax/5);
    flux_q3(t)     = Q3(flux_q2(t),tmp_sFlux(3),smax/5);
    flux_q4(t)     = Q4(flux_q3(t),tmp_sFlux(4),smax/5);
    flux_r3(t)     = R3(flux_q4(t),tmp_sFlux(5),smax/5);
    flux_rg(t)     = RG(g,flux_r3(t));
    flux_r3star(t) = R3STAR(1-g,flux_r3(t));
    flux_qg(t)     = QG(kg,tmp_sFlux(6));
        
    % Update the stores
    store_S1(t) = S1old + (flux_i(t)  - flux_e1(t) - flux_q1(t)) * delta_t;
    store_S2(t) = S2old + (flux_q1(t) - flux_e2(t) - flux_q2(t)) * delta_t;    
    store_S3(t) = S3old + (flux_q2(t) - flux_e3(t) - flux_q3(t)) * delta_t;
    store_S4(t) = S4old + (flux_q3(t) - flux_e4(t) - flux_q4(t)) * delta_t;
    store_S5(t) = S5old + (flux_q4(t) - flux_e5(t) - flux_r3(t)) * delta_t;
    store_S6(t) = S6old + (flux_rg(t) - flux_qg(t)) * delta_t; 

% Routing -----------------------------------------------------------------
    tmp_Qr_cur    = (flux_r1(t) + flux_r2(t) + flux_r3star(t)).*uh_gamma;   % find how the runoff of this time step will be distributed in time
    tmp_Qr_old    = tmp_Qr_old + tmp_Qr_cur;                                % update the 'still-to-flow-out' vector
    flux_qr(t)    = tmp_Qr_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_Qr_old    = circshift(tmp_Qr_old,-1);                               % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_Qr_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)
    
    % The Gamma-based routing function technically extends from 0 to
    % infinity. The UH is only calcualted for a duration equal to the
    % maximum simulation length. This variable tracks water routed beyond
    % this time, which we need for water balance calculation.
    tmp_Qr_beyond_tEnd(t) = frac_routing_beyond_time_series * (flux_r1(t) + flux_r2(t) + flux_r3star(t));

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_evap + flux_e1 + flux_e2 + flux_e3 + flux_e4 + flux_e5) * delta_t;
    fluxOutput.Q      = (flux_qr + flux_qg) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pstar  = flux_pstar * delta_t;
    fluxInternal.r1     = flux_r1 * delta_t;
    fluxInternal.i      = flux_i * delta_t;
    fluxInternal.r2     = flux_r2 * delta_t;
    fluxInternal.e1     = flux_e1 * delta_t;
    fluxInternal.e2     = flux_e2 * delta_t;
    fluxInternal.e3     = flux_e3 * delta_t;
    fluxInternal.e4     = flux_e4 * delta_t;
    fluxInternal.e5     = flux_e5 * delta_t;
    fluxInternal.q1     = flux_q1 * delta_t;
    fluxInternal.q2     = flux_q2 * delta_t;
    fluxInternal.q3     = flux_q3 * delta_t;
    fluxInternal.q4     = flux_q4 * delta_t;
    fluxInternal.r3     = flux_r3 * delta_t;
    fluxInternal.rg     = flux_rg * delta_t;
    fluxInternal.r3star = flux_r3star * delta_t;
    fluxInternal.qg     = flux_qg * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;

% Check water balance
if nargout == 4
    
    % Add the 'routing-beyond-simulation-time' variable to the
    % 'still-to-be-routed' variable
    tmp_rout = tmp_Qr_beyond_tEnd;
    tmp_rout(1:length(tmp_Qr_old)) = tmp_rout(1:length(tmp_Qr_old)) + tmp_Qr_old;
    
    % Calculate water balance
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     tmp_rout);         % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end

 




