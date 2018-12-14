function [ fluxOutput, fluxInternal, storeInternal, waterBalance  ] = ...
            m_13_hillslope_7p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Hillslope model 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Savenije, H. H. G. (2010). “Topography driven conceptual modelling 
% (FLEX-Topo).” Hydrology and Earth System Sciences, 14(12), 2681–2692. 
% https://doi.org/10.5194/hess-14-2681-2010


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
t_end = length(P);

% Parameters 
% [name in documentation] = theta(order in which specified in parameter file)
dw    = theta(1);     % Daily interception [mm]
betaw = theta(2);     % Soil moisture storage distribution parameter [-]
swmax = theta(3);     % Maximum soil moisture storage [mm]
a     = theta(4);     % Division parameter for surface and groundwater flow [-]
th    = theta(5);     % Time delay of surface flow [d]
c     = theta(6);     % Rate of capillary rise [mm/d]
kh    = theta(7);     % Groundwater runoff coefficient [d-1]

%%INITIALISE MODEL STORES
S10         = storeInitial(1);      % Initial soil moisture storage
S20         = storeInitial(2);      % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                  % lower bounds of stores
store_upp = [];                     % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_pe    = zeros(1,t_end);
flux_ei    = zeros(1,t_end);
flux_ea    = zeros(1,t_end);
flux_qse   = zeros(1,t_end);
flux_qses  = zeros(1,t_end);
flux_qseg  = zeros(1,t_end);
flux_c     = zeros(1,t_end);
flux_qhsrf = zeros(1,t_end);
flux_qhgw  = zeros(1,t_end);
flux_qt    = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_half] = uh_3_half(1,th,delta_t);

%%INITIALISE ROUTING VECTOR
Qhsrf_old   = zeros(1,length(uh_half));    % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture
% S2. Groundwater

% PE(P(t),dw): effective rainfall after interception
PE = interception_2;

% EA(S1,Ep(t),delta_t): evaporation from soil moisture
EA = evap_1;

% QSE(S1,swmax,betaw,PE(P(t),dw)): saturation excess flow
QSE = saturation_2;

% C(c,S2,delta_t): capillary rise
C = capillary_2;

% QSES(a,QSE(S1,swmax,betaw,PE)): fraction saturation excess to stream
QSES = split_1;

% QSEG(1-a,QSE(S1,swmax,betaw,PE)): fraction saturation excess to groundwater
QSEG = split_1;

% QHGW(kh,S2): baseflow
QHGW = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % Disable display settings
                                'JacobPattern', [1,1;
                                                 1,1]);                     % Specify the Jacobian pattern                                             
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
    tmpf_S1 = @(S1,S2) (PE(P(t),dw) + ...
                        C(c,S2,delta_t) - ...
                        EA(S1,Ep(t),delta_t) - ...
                        QSE(S1,swmax,betaw,PE(P(t),dw)));                   % store 1 function
    tmpf_S2 = @(S1,S2) ((1-a).*QSE(S1,swmax,betaw,PE(P(t),dw)) - ...
                        C(c,S2,delta_t) - ...
                        QHGW(kh,S2));                                       % store 2 function
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2)),...                            % system of storage equations
                        [S1old,S2old],...                                   % storage values on previous time step
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
                                                    eq_sys(1),eq_sys(2)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old], ...                  % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
                    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_pe(t)   = PE(P(t),dw);
    flux_ei(t)   = P(t)-flux_pe(t);                                         % tracks 'intercepted' rainfall
    flux_ea(t)   = EA(tmp_sFlux(1),Ep(t),delta_t);
    flux_qse(t)  = QSE(tmp_sFlux(1),swmax,betaw,flux_pe(t));
    flux_qses(t) = QSES(a,flux_qse(t));
    flux_qseg(t) = QSEG(1-a,flux_qse(t));
    flux_c(t)    = C(c,tmp_sFlux(2),delta_t);
    flux_qhgw(t) = QHGW(kh,tmp_sFlux(2));

    % Do the surface routing
    Qhrs_cur      = flux_qses(t).*uh_half;                                  % find how the runoff of this time step will be distributed in time
    Qhsrf_old     = Qhsrf_old + Qhrs_cur;                                   % update the 'still-to-flow-out' vector
    flux_qhsrf(t) = Qhsrf_old(1);                                           % the first value in 'still-to-flow-out' vector leaves the model this time step
    Qhsrf_old = circshift(Qhsrf_old,-1);                                    % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    Qhsrf_old(end) = 0;                                                     % update the 'still-to-flow-out' vector
    
    % total outflow
    flux_qt(t)  = flux_qhsrf(t) + flux_qhgw(t);
    
    % Update the stores
    store_S1(t) = S1old + (flux_pe(t) + flux_c(t) - flux_ea(t) - flux_qse(t)) * delta_t;
    store_S2(t) = S2old + (flux_qseg(t) - flux_c(t) - flux_qhgw(t))  * delta_t;

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = (flux_ea + flux_ei) * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pe    = flux_pe * delta_t;
    fluxInternal.ei    = flux_ei;
    fluxInternal.ea    = flux_ea;
    fluxInternal.qse   = flux_qse * delta_t;
    fluxInternal.qseg  = flux_qseg * delta_t;
    fluxInternal.qses  = flux_qses * delta_t;
    fluxInternal.c     = flux_c * delta_t;
    fluxInternal.Qhsrf = flux_qhsrf * delta_t;
    fluxInternal.Qhgw  = flux_qhgw * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     Qhsrf_old);        % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end



 




