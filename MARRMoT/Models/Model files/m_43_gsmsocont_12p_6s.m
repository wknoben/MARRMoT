function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
        m_43_gsmsocont_12p_6s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [GSM-SOCONT] 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Schaefli, B., Hingray, B., Niggli, M., & Musy, a. (2005). A conceptual 
% glacio-hydrological model for high mountainous catchments. Hydrology and 
% Earth System Sciences Discussions, 2(1), 73–117. 
% http://doi.org/10.5194/hessd-2-73-2005


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
fice  = theta(1);     % Fraction of catchment covered by glacier [-]
t0    = theta(2);     % Threshold temperature for snowfall [oC]
asnow = theta(3);     % Degree-day factor for snow melt [mm/oC/d]
tm    = theta(4);     % Threshold temperature for snow melt [oC]
ks    = theta(5);     % Runoff coeficient for snow melt on glacier [d-1]
aice  = theta(6);     % Threshold temperature for ice melt [oC]
ki    = theta(7);     % Runoff coeficient for ice melt on glacier [d-1]
a     = theta(8);     % Maximum soil moisture storage [mm]
x     = theta(9);     % Evaporation non-linearity [-]
y     = theta(10);    % Infiltration non-linearity [-]
ksl   = theta(11);    % Runoff coefficient for baseflow [d-1]
beta  = theta(12);    % Runoff coefficient for quick flow [mm^(4/3)/d]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial snow on glacier storage
S20   = storeInitial(2);       % Initial snow on glacier runoff storage
S30   = storeInitial(3);       % Initial ice on glacier runoff storage
S40   = storeInitial(4);       % Initial snow on soil storage
S50   = storeInitial(5);       % Initial soil moisture storage
S60   = storeInitial(6);       % Initial quick flow storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0];     % lower bounds of stores
store_upp = [];                % optional higher bounds

%%INITIALISE STORAGE VECTORS 
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);

flux_pice  = zeros(1,t_end);
flux_pices = zeros(1,t_end);
flux_picer = zeros(1,t_end);
flux_mis   = zeros(1,t_end);
flux_pirs  = zeros(1,t_end);
flux_piri  = zeros(1,t_end);
flux_mii   = zeros(1,t_end);
flux_qis   = zeros(1,t_end);
flux_qii   = zeros(1,t_end);
flux_pni   = zeros(1,t_end);
flux_pnis  = zeros(1,t_end);
flux_pnir  = zeros(1,t_end);
flux_mnis  = zeros(1,t_end);
flux_peq   = zeros(1,t_end);
flux_peff  = zeros(1,t_end);
flux_pinf  = zeros(1,t_end);
flux_et    = zeros(1,t_end);
flux_qsl   = zeros(1,t_end);
flux_qqu   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Snow on glacier storage
% S2. Snow on glacier runoff storage
% S3. Ice on glacier runoff storage
% S4. Snow on soil storage
% S5. Soil moisture storage
% S6. Quick flow storage

% PICE(fice,P(t)): fraction of P on glacier
PICE = split_1;

% PICES(PICE(fice,P(t)),T(t),t0): snowfall on glacier
PICES = snowfall_1;

% PICER(PICE(fice,P(t)),T(t),t0): rainfall on glacier
PICER = rainfall_1;

% MIS(asnow,tm,T(t),S1,delta_t): snowmelt on glacier
MIS = melt_1;

% PIRS(PICER(PICE(fice,P(t)),T(t),t0),S1,0.01): rainfall on snow 
PIRS = saturation_9;

% PIRI(PICER(PICE(fice,P(t)),T(t),t0),PIRS(PICER(PICE(fice,P(t)),T(t),t0),S1,0.01)): 
% rainfall on ice
PIRI = effective_1;

% MII(aice,tm,T(t),Inf,S1,0.01,delta_t): glacier melt conditional on no snow pack
MII = melt_3;

% QIS(ks,S2): flow from snow melt on glacier
QIS = baseflow_1;

% QII(ki,S3): flow from ice melt on glacier
QII = baseflow_1;

% PNI(1-fice,P(t)): fraction of P on soil
PNI = split_1;

% PNIS(PNI(1-fice,P(t)),T(t),t0): snowfall on soil
PNIS = snowfall_1;

% PNIR(PNI(1-fice,P(t)),T(t),t0): rainfall on soil
PNIR = rainfall_1;

% MNIS(asnow,tm,T(t),S4,delta_t): snowmelt on soil
MNIS = melt_1;

% PEFF(1,y,S5,a,PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4)):
% infiltration into soil
PEFF = infiltration_6;

% PINF(PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4),PEFF(1,y,S5,a,PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4))):
% direct runoff
PINF = effective_1;

% ET(1,x,S5,a,Ep(t),delta_t): evapotranspiration
ET = evap_19;

% QSL(ksl,S5): slow flow
QSL = baseflow_1;

% QQU(beta,5/3,S6,delta_t): quick flow
QQU = interflow_3;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0;
                                               1,1,0,0,0,0;
                                               1,0,1,0,0,0;
                                               0,0,0,1,0,0;
                                               0,0,0,1,1,0;
                                               0,0,0,1,1,1]);               % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,0,0,0,0,0;
                                                    1,1,0,0,0,0;
                                                    1,0,1,0,0,0;
                                                    0,0,0,1,0,0;
                                                    0,0,0,1,1,0;
                                                    0,0,0,1,1,1],...
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
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5,S6) ...
               (PICES(PICE(fice,P(t)),T(t),t0) - ...
                MIS(asnow,tm,T(t),S1,delta_t));
    
    tmpf_S2 = @(S1,S2,S3,S4,S5,S6) ...
               (PIRS(PICER(PICE(fice,P(t)),T(t),t0),S1,0.01) + ...
                MIS(asnow,tm,T(t),S1,delta_t) - ...
                QIS(ks,S2));
           
    tmpf_S3 = @(S1,S2,S3,S4,S5,S6) ...
               (PIRI(PICER(PICE(fice,P(t)),T(t),t0),PIRS(PICER(PICE(fice,P(t)),T(t),t0),S1,0.01)) + ...
                MII(aice,tm,T(t),Inf,S1,0.01,delta_t) - ...
                QII(ki,S3));
         
    tmpf_S4 = @(S1,S2,S3,S4,S5,S6) ...
               (PNIS(PNI(1-fice,P(t)),T(t),t0) - ...
                MNIS(asnow,tm,T(t),S4,delta_t));
           
    tmpf_S5 = @(S1,S2,S3,S4,S5,S6) ...
               (PINF(PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4,delta_t),PEFF(1,y,S5,a,PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4,delta_t))) - ...
                ET(1,x,S5,a,Ep(t),delta_t) - ...
                QSL(ksl,S5));
           
	tmpf_S6 = @(S1,S2,S3,S4,S5,S6) ...
               (PEFF(1,y,S5,a,PNIR(PNI(1-fice,P(t)),T(t),t0)+MNIS(asnow,tm,T(t),S4,delta_t)) - ...
                QQU(beta,5/3,S6,delta_t));
           
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old,S6old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5,tmpf_S6);     % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),...
                        eq_sys(4),eq_sys(5),eq_sys(6)),...                  % system of storage equations
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
                                        eq_sys(1),eq_sys(2),eq_sys(3),...
                                        eq_sys(4),eq_sys(5),eq_sys(6)),...
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
    flux_pice(t)    = PICE(fice,P(t));
    flux_pices(t)   = PICES(flux_pice(t),T(t),t0);
    flux_picer(t)   = PICER(flux_pice(t),T(t),t0);
    flux_mis(t)     = MIS(asnow,tm,T(t),tmp_sFlux(1),delta_t);
    flux_pirs(t)    = PIRS(flux_picer(t),tmp_sFlux(1),0.01);
    flux_piri(t)    = PIRI(flux_picer(t),flux_pirs(t));
    flux_mii(t)     = MII(aice,tm,T(t),Inf,tmp_sFlux(1),0.01,delta_t);
    flux_qis(t)     = QIS(ks,tmp_sFlux(2));
    flux_qii(t)     = QII(ki,tmp_sFlux(3));
    flux_pni(t)     = PNI(1-fice,P(t));
    flux_pnis(t)    = PNIS(flux_pni(t),T(t),t0);
    flux_pnir(t)    = PNIR(flux_pni(t),T(t),t0);
    flux_mnis(t)    = MNIS(asnow,tm,T(t),tmp_sFlux(4),delta_t);
    flux_peq(t)     = flux_pnir(t) + flux_mnis(t);
    flux_peff(t)    = PEFF(1,y,tmp_sFlux(5),a,flux_peq(t));
    flux_pinf(t)    = PINF(flux_peq(t),flux_peff(t));
    flux_et(t)      = ET(1,x,tmp_sFlux(5),a,Ep(t),delta_t);
    flux_qsl(t)     = QSL(ksl,tmp_sFlux(5));
    flux_qqu(t)     = QQU(beta,5/3,tmp_sFlux(6),delta_t);
    
    % Update the stores
    store_S1(t) = S1old + (flux_pices(t) - flux_mis(t)) * delta_t;
    store_S2(t) = S2old + (flux_pirs(t)  + flux_mis(t) - flux_qis(t)) * delta_t;
    store_S3(t) = S3old + (flux_piri(t)  + flux_mii(t) - flux_qii(t)) * delta_t;
    store_S4(t) = S4old + (flux_pnis(t)  - flux_mnis(t)) * delta_t;
    store_S5(t) = S5old + (flux_pinf(t)  - flux_et(t) - flux_qsl(t)) * delta_t;  
    store_S6(t) = S6old + (flux_peff(t)  - flux_qqu(t)) * delta_t;
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_et) * delta_t;
    fluxOutput.Q      = (flux_qii + flux_qis + flux_qsl + flux_qqu) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pice  = flux_pice * delta_t;
    fluxInternal.pices = flux_pices * delta_t;
    fluxInternal.picer = flux_picer * delta_t;
    fluxInternal.mis   = flux_mis * delta_t;
    fluxInternal.pirs  = flux_pirs * delta_t;
    fluxInternal.piri  = flux_piri * delta_t;
    fluxInternal.mii   = flux_mii * delta_t;
    fluxInternal.qis   = flux_qis * delta_t;
    fluxInternal.qii   = flux_qii * delta_t;
    fluxInternal.pni   = flux_pni * delta_t;
    fluxInternal.pnis  = flux_pnis * delta_t;
    fluxInternal.pnir  = flux_pnir * delta_t;
    fluxInternal.mnis  = flux_mnis * delta_t;
    fluxInternal.peq   = flux_peq * delta_t;
    fluxInternal.peff  = flux_peff * delta_t;
    fluxInternal.pinf  = flux_pinf * delta_t;
    fluxInternal.et    = flux_et * delta_t;
    fluxInternal.qsl   = flux_qsl * delta_t;
    fluxInternal.qqu   = flux_qqu * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;

% Check water balance
if nargout == 4
    
    % Manual water balance on account of the infinite glacier depth
    waterBalance =  sum(P).*delta_t + ...       % Incoming precipitation
                    sum(flux_mii).*delta_t - ...% Incoming glacier melt
                    sum(fluxOutput.Q) - ...     % Outgoing flow
                    sum(fluxOutput.Ea) - ...    % Outgoing evaporation
                    (store_S1(end)-S10) - ...   % Storage change
                    (store_S2(end)-S20) - ...
                    (store_S3(end)-S30) - ...
                    (store_S4(end)-S40) - ...
                    (store_S5(end)-S50) - ...
                    (store_S6(end)-S60);
    
    % Add the glacier melt flux to the water balance
    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Glacier melt = ',num2str(sum(flux_mii).*delta_t),' mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['Water balance = sum(P) + sum(glacier melt) - (sum(Q) + sum(E_a)) - delta S = ',num2str(waterBalance)]);  

end



 




