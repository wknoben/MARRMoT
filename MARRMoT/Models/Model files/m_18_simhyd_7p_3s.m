function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
    m_18_simhyd_7p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: SIMHYD 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Chiew, F. H. S., Peel, M. C., & Western, A. W. (2002). Application and 
% testing of the simple rainfall-runoff model SIMHYD. In V. P. Singh & D. 
% K. Frevert (Eds.), Mathematical Models of Small Watershed Hydrology (pp. 
% 335–367). Chelsea, Michigan, USA: Water Resources Publications LLC, USA.


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
insc    = theta(1);     % Maximum interception capacity, [mm]
coeff   = theta(2);     % Maximum infiltration loss parameter, [-]
sq      = theta(3);     % Infiltration loss exponent, [-]
smsc    = theta(4);     % Maximum soil moisture capacity, [mm]
sub     = theta(5);     % Proportionality constant, [-],
crak    = theta(6);     % Proportionality constant, [-]
k       = theta(7);     % Slow flow time scale, [d-1]

%%INITIALISE MODEL STORES 
S10 = storeInitial(1);       % Initial interception storage
S20 = storeInitial(2);       % Initial soil moisture storage
S30 = storeInitial(3);       % Initial groundwater storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];         % lower bounds of stores
store_upp = [];              % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1     = zeros(1,t_end);     
store_S2     = zeros(1,t_end);     
store_S3     = zeros(1,t_end);     
flux_Ei      = zeros(1,t_end);     
flux_EXC     = zeros(1,t_end);     
flux_INF     = zeros(1,t_end);     
flux_INT     = zeros(1,t_end); 
flux_REC     = zeros(1,t_end); 
flux_Et      = zeros(1,t_end);
flux_GWF     = zeros(1,t_end);
flux_BAS     = zeros(1,t_end);
flux_SRUN    = zeros(1,t_end);
flux_Qt      = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Soil moisture
% S3. Groundwater

% Ei(S1,Ep(t),delta_t): evaporation from interception
Ei = evap_1;    

% EXC(P(t),S1,insc): overflow from interception
EXC = interception_1; 

% INF(coeff,sq,S2,smsc,EXC): infiltration of interception excess
INF = infiltration_1;

% INT(sub,S2,smsc,INF): interflow and saturation excess flow
INT = interflow_1;

% REC(crak,S2,smsc,INF-INT): recharge to groundwater
REC = recharge_1;

% Et(10,S2,smsc,Ep(t),delta_t): scaled evaporation from soil moisture
Et = evap_2;

% GWF(INF-INT-REC,S2,smsc): overflow from soil moisture
GWF = saturation_1;

% BAS(k,S3): baseflow from groundwater. 
BAS = baseflow_1;

% SRUN(EXC,INF): surface runoff after infiltration.
%SRUN is defined when the final fluxes are calculated

% Qt(SRUN,INT,BAS): total outflow from surface, interflow and baseflow.
% Qt is defined when the final fluxes are calculated

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % Disable display settings
                                'JacobPattern', [1,0,0;
                                                 1,1,0;
                                                 1,1,1]);                   % Specify the Jacobian pattern
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases where fsolve fails
                                 'Display','none',...
                                 'JacobPattern',[1,0,0;
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
        % Interception store
        tmpf_S1 = @(S1,S2,S3) (P(t) - ...
                               Ei(S1,Ep(t),delta_t) - ...
                               EXC(P(t),S1,insc));
    
        % Soil moisture
        tmpf_S2 = @(S1,S2,S3) (INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                                INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                                REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))) - ...
                               Et(10,S2,smsc,Ep(t),delta_t) - ...
                               GWF((INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                                INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                                REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))),...
                                S2,smsc);
                               
        % Groundwater
        tmpf_S3 = @(S1,S2,S3) REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))))) + ...
                              GWF((INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                                INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                                REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))),...
                                S2,smsc) - ...
                              BAS(k,S3);
        
    
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
    % Calculate the fluxes
    flux_Ei(t)   = Ei(tmp_sNew(1),Ep(t),delta_t);
    flux_EXC(t)  = EXC(P(t),tmp_sNew(1),insc);
    flux_INF(t)  = INF(coeff,sq,tmp_sNew(2),smsc,flux_EXC(t));
    flux_INT(t)  = INT(sub,tmp_sNew(2),smsc,flux_INF(t));
    flux_REC(t)  = REC(crak,tmp_sNew(2),smsc,(flux_INF(t)-flux_INT(t)));
    flux_Et(t)   = Et(10,tmp_sNew(2),smsc,Ep(t),delta_t);
    flux_GWF(t)  = GWF((flux_INF(t)-flux_INT(t)-flux_REC(t)),tmp_sNew(2),smsc);
    flux_BAS(t)  = BAS(k,tmp_sNew(3));
    flux_SRUN(t) = flux_EXC(t) - flux_INF(t); 
    flux_Qt(t)   = flux_SRUN(t) + flux_INT(t) + flux_BAS(t);
    
    % update the stores
    store_S1(t)  = S1old + (P(t) - flux_Ei(t) - flux_EXC(t)) * delta_t;
    store_S2(t)  = S2old + ((flux_INF(t) - flux_INT(t) - flux_REC(t)) - ...
                            flux_Et(t) - flux_GWF(t)) * delta_t;   
    store_S3(t)  = S3old + (flux_REC(t) + flux_GWF(t) - flux_BAS(t)) * delta_t;
    
end

%% 6. Generate outputs
% --- Fluxes leaving the model ---
    fluxOutput.Ea     = (flux_Ei + flux_Et) * delta_t;
    fluxOutput.Q      = flux_Qt * delta_t;
        
    % --- Fluxes internal to the model ---
    fluxInternal.EXC  = flux_EXC * delta_t;
    fluxInternal.INF  = flux_INF * delta_t;
    fluxInternal.INT  = flux_INT * delta_t;
    fluxInternal.REC  = flux_REC * delta_t;
    fluxInternal.GWF  = flux_GWF * delta_t;
    fluxInternal.BAS  = flux_BAS * delta_t;
    fluxInternal.SRUN = flux_SRUN * delta_t;
    
    % --- Stores ---
    storeInternal.S1 = store_S1;
    storeInternal.S2 = store_S2;
    storeInternal.S3 = store_S3;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end   



 




