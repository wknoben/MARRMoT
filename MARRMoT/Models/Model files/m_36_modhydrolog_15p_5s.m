function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
        m_36_modhydrolog_15p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: MODHYDROLOG 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Chiew, F. H. S. (1990). Estimating groundwater recharge using an
% integrated surface and groundwater model. University of Melbourne.
%
% Chiew, F., & McMahon, T. (1994). Application of the daily rainfall-runoff
% model MODHYDROLOG to 28 Australian catchments. Journal of Hydrology, 
% 153(1–4), 383–416. https://doi.org/10.1016/0022-1694(94)90200-3

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
sub     = theta(5);     % Proportionality constant, [-]
crak    = theta(6);     % Proportionality constant, [-]
em      = theta(7);     % Plant-controled maximum evaporation rate [mm/d]
dsc     = theta(8);     % Maximum depression capacity, [mm]
ads     = theta(9);     % Land fraction functioning as depression storage, [-]
md      = theta(10);    % Depression storage parameter, [-], default = 1
vcond   = theta(11);    % Leakage coefficient, [mm/d]
dlev    = theta(12);    % Datum around which groundwater fluctuates relative to river bed, [mm]
k1      = theta(13);    % Flow exchange parameter, [d-1] 
k2      = theta(14);    % Flow exchange parameter, [d-1] 
k3      = theta(15);    % Flow exchange parameter, [d-1] 

%%INITIALISE MODEL STORES AND ROUTING VECTORS
S10 = storeInitial(1);  % Initial interception storage
S20 = storeInitial(2);  % Initial soil moisture storage
S30 = storeInitial(3);  % Initial depression storage
S40 = storeInitial(4);  % Initial groundwater storage
S50 = storeInitial(5);  % Initial channel storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,-inf,0];  % lower bounds of stores
store_upp = [];              % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);

flux_Ei     = zeros(1,t_end);         
flux_EXC    = zeros(1,t_end);
flux_INF    = zeros(1,t_end);
flux_INT    = zeros(1,t_end);
flux_REC    = zeros(1,t_end);
flux_SMF    = zeros(1,t_end);
flux_Et     = zeros(1,t_end);
flux_GWF    = zeros(1,t_end);
flux_TRAP   = zeros(1,t_end);
flux_Ed     = zeros(1,t_end);
flux_DINF   = zeros(1,t_end);
flux_SEEP   = zeros(1,t_end);
flux_FLOW   = zeros(1,t_end);
flux_Q      = zeros(1,t_end);
flux_RUN    = zeros(1,t_end);
flux_SRUN   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Soil moisture
% S3. Depression
% S4. Groundwater
% S5. Channel

% EI(S1,Ep,delta_t): evaporation from interception
EI = evap_1;    

% EXC(P,S1,insc): overflow from interception
EXC = interception_1; 

% INF(coeff,sq,S2,smsc,EXC): infiltration of interception excess
INF = infiltration_1;

% INT(sub,S2,smsc,INF): interflow and saturation excess flow
INT = interflow_1;

% REC(crak,S2,smsc,INF-INT): recharge to groundwater
REC = recharge_1;

% ET(em,S2,smsc,Ep,delta_t): scaled evaporation from soil moisture
ET = evap_2;

% GWF(INF-INT-REC,S2,smsc): overflow from soil moisture
GWF = saturation_1; 

% TRAP(ads,md,S3,dsc,EXC-INF,delta_t): storage of surface flow in depressions
TRAP = depression_1;

% ED(S3,ads*Ep,delta_t): evaporation from depression
ED = evap_1;

% DINF(ads,coeff,sq,S2,smsc,[INF-INT-REC],S3)
tmp_dinf = infiltration_2;                                                  % needed because we need to aply a scaling coefficient to this equation
DINF = @(ads,coeff,sq,S2,smsc,IN,S3,delta_t) ads.*tmp_dinf(coeff,sq,S2,smsc,IN,S3,delta_t);
clear tmp_dinf                                                              % no longer needed

% SEEP(vcond,S4,dlev): seepage to or from deeper aquifer 
SEEP = exchange_3;

% FLOW(k1,k2,k3,S4,SRUN,delta_t): exchange between river and groundwater
FLOW = exchange_1;

% Q(1,S5): outflow from channel
Q = baseflow_1;                                                             % special case of the linear reservoir, Q = 1*S

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...
                                'JacobPattern', [1,0,0,0,0;
                                                 1,1,0,0,0;
                                                 1,1,1,0,0;
                                                 1,1,1,1,0;
                                                 1,1,1,1,1]);
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0,0;
                                                  1,1,0,0,0;
                                                  1,1,1,0,0;
                                                  1,1,1,1,0;
                                                  1,1,1,1,1],...
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
                    (P(t) - ...
                     EI(S1,Ep(t),delta_t) - ...
                     EXC(P(t),S1,insc));
                           
        tmpf_S2 = @(S1,S2,S3,S4,S5) ...
                    (INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                        INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                        REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))) + ...
                     DINF(ads,coeff,sq,S2,smsc, ...
                        (INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                        INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                        REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))), ...
                        S3,delta_t) - ...   
                     ET(em,S2,smsc,Ep(t),delta_t) - ...
                     GWF((INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                        INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                        REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))),...
                        S2,smsc);
        
        tmpf_S3 = @(S1,S2,S3,S4,S5) ...
                    TRAP(ads,md,S3,dsc,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))),delta_t) - ...
                    ED(S3,ads*Ep(t),delta_t) - ...
                    DINF(ads,coeff,sq,S2,smsc,...
                        (INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                        INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                        REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))),...
                        S3,delta_t);
        
        tmpf_S4 = @(S1,S2,S3,S4,S5) ...
                    REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))))) + ...
                    GWF((INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)) - ...
                        INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...
                        REC(crak,S2,smsc,(INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))-INT(sub,S2,smsc,INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc)))))),...
                        S2,smsc) - ...
                    SEEP(vcond,S4,dlev) - ...
                    FLOW(k1,k2,k3,S4,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...        
                        TRAP(ads,md,S3,dsc,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))),delta_t),delta_t);
                    
        tmpf_S5 = @(S1,S2,S3,S4,S5) ...
                    (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...        
                        TRAP(ads,md,S3,dsc,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))),delta_t) + ...
                    INT(sub,S2,smsc,...
                        INF(coeff,sq,S2,smsc,...
                        EXC(P(t),S1,insc))) + ...
                    FLOW(k1,k2,k3,S4,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))) - ...        
                        TRAP(ads,md,S3,dsc,...
                        (EXC(P(t),S1,insc) - INF(coeff,sq,S2,smsc,EXC(P(t),S1,insc))),delta_t),delta_t) - ...
                    Q(1,S5);
                
    % Call the numerical scheme function to create the ODE approximations.
    % This returns a new anonymous function that we solve in the next step
    solve_fun = feval(scheme,...                                            % numerical scheme
                    [S1old,S2old,S3old,S4old,S5old],...                     % old storages
                    delta_t,...                                             % time step
                    tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5);               % temporary ODE's
    
% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) ...
                        solve_fun(eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                            eq_sys(5)),...                                  % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old],...                 % storage values on previous time step
                        fsolve_options);            
        
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
                                        eq_sys(4),eq_sys(5)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old,S5old], ...% storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds
                                              
    end

%% End of Backup of simultaneous solving
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_Ei(t)     = EI(tmp_sFlux(1),Ep(t),delta_t);         
    flux_EXC(t)    = EXC(P(t),tmp_sFlux(1),insc);
    flux_INF(t)    = INF(coeff,sq,tmp_sFlux(2),smsc,flux_EXC(t));
    flux_INT(t)    = INT(sub,tmp_sFlux(2),smsc,flux_INF(t));
    flux_REC(t)    = REC(crak,tmp_sFlux(2),smsc,flux_INF(t)-flux_INT(t));
    flux_SMF(t)    = flux_INF(t) - flux_INT(t) - flux_REC(t);
    flux_Et(t)     = ET(em,tmp_sFlux(2),smsc,Ep(t),delta_t);
    flux_GWF(t)    = GWF(flux_SMF(t),tmp_sFlux(2),smsc);
    flux_RUN(t)    = flux_EXC(t) - flux_INF(t);
    flux_TRAP(t)   = TRAP(ads,md,tmp_sFlux(3),dsc,flux_RUN(t),delta_t);
    flux_Ed(t)     = ED(tmp_sFlux(3),ads*Ep(t),delta_t);
    flux_DINF(t)   = DINF(ads,coeff,sq,tmp_sFlux(2),smsc,flux_SMF(t),tmp_sFlux(3),delta_t);
    flux_SEEP(t)   = SEEP(vcond,tmp_sFlux(4),dlev);
    flux_SRUN(t)   = flux_RUN(t) - flux_TRAP(t);
    flux_FLOW(t)   = FLOW(k1,k2,k3,tmp_sFlux(4),flux_SRUN(t),delta_t);
    flux_Q(t)      = Q(1,tmp_sFlux(5));
    
    % update the stores
    store_S1(t) = S1old + (P(t)          - flux_Ei(t)    - flux_EXC(t)) * delta_t;
    store_S2(t) = S2old + (flux_SMF(t)   + flux_DINF(t)  - flux_Et(t)    - flux_GWF(t)) * delta_t;
    store_S3(t) = S3old + (flux_TRAP(t)  - flux_Ed(t)    - flux_DINF(t)) * delta_t;
    store_S4(t) = S4old + (flux_REC(t)   + flux_GWF(t)   - flux_SEEP(t)  - flux_FLOW(t)) * delta_t;
    store_S5(t) = S5old + (flux_SRUN(t)  + flux_INT(t)   + flux_FLOW(t)  - flux_Q(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = (flux_Ei + flux_Et + flux_Ed) * delta_t;
    fluxOutput.Q      = flux_Q * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.Ei   = flux_Ei * delta_t;
    fluxInternal.Et   = flux_Et * delta_t;
    fluxInternal.Ed   = flux_Ed * delta_t; 
    fluxInternal.exc  = flux_EXC * delta_t;
    fluxInternal.inf  = flux_INF * delta_t;
    fluxInternal.int  = flux_INT * delta_t;
    fluxInternal.rec  = flux_REC * delta_t;
    fluxInternal.smf  = flux_SMF * delta_t;
    fluxInternal.run  = flux_RUN * delta_t;
    fluxInternal.trap = flux_TRAP * delta_t;
    fluxInternal.srun = flux_SRUN * delta_t;
    fluxInternal.dinf = flux_DINF * delta_t;
    fluxInternal.gwf  = flux_GWF * delta_t;
    fluxInternal.seep = flux_SEEP * delta_t;
    fluxInternal.flow = flux_FLOW * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;

% Check water balance
if nargout == 4
    waterBalance =  sum(P).*delta_t - ...               % Incoming precipitation
                    sum(fluxOutput.Q) - ...             % Outgoing flow
                    sum(fluxOutput.Ea) - ...            % Outgoing evaporation
                    (store_S1(end)-S10) - ...           % Storage change
                    (store_S2(end)-S20) - ...
                    (store_S3(end)-S30) - ...
                    (store_S4(end)-S40) - ...
                    (store_S5(end)-S50) - ...
                    sum(fluxInternal.seep);             % Seepage to deeper groundwater

    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
    disp(['Delta S4 = ',num2str((store_S4(end)-S40)),'mm.'])
    disp(['Delta S5 = ',num2str((store_S5(end)-S50)),'mm.'])
    disp(['Water balance = sum(P) - sum(SEEP) - (sum(Q) + sum(E_a)) - delta S = ',num2str(waterBalance)]);

end


