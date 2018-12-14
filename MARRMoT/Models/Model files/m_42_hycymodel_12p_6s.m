function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_42_hycymodel_12p_6s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: HYCYMODEL 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Fukushima, Y. (1988). A model of river flow forecasting for a small 
% forested mountain catchment. Hydrological Processes, 2(2), 167–185.

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
c       = theta(1);     % Fraction area that is channel [-]
imax    = theta(2);     % Maximum total interception storage [mm]
a       = theta(3);     % Fraction stem/trunk interception [-]
fi2     = theta(4);     % Fraction of total interception that is trunk/stem interception [mm]
kin     = theta(5);     % Infiltration runoff coefficient [d-1]
d50     = theta(6);     % Soil depth where 50% of area contributes to effective flow [mm]
fd16    = theta(7);     % Fraction of D50 that is D16 [-]
sbc     = theta(8);     % Soil depth where evaporation rate starts to decline [mm]
kb      = theta(9);     % Baseflow runoff coefficient [d-1]
pb      = theta(10);    % Baseflow non-linearity [-]
kh      = theta(11);    % Hillslope runoff coefficient [d-1]
kc      = theta(12);    % Channel runoff coefficient [d-1]

% Auxiliary parameters 
i1max   = (1-fi2)*imax; % Maximum canopy interception [mm]
i2max   = fi2*imax;     % Maximum trunk/stem interception [mm]
d16     = fd16*d50;     % Soil depth where 16% of area contributes to effective flow [mm]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial canopy storage
S20     = storeInitial(2);       % Initial stem/trunk storage
S30     = storeInitial(3);       % Initial upper soil moisture storage
S40     = storeInitial(4);       % Initial baseflow routing
S50     = storeInitial(5);       % Initial hillslope routing
S60     = storeInitial(6);       % Initial channel routing

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0];       % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);

flux_rc  = zeros(1,t_end);
flux_rg  = zeros(1,t_end);
flux_eic = zeros(1,t_end);
flux_qie = zeros(1,t_end);
flux_qis = zeros(1,t_end);
flux_rt  = zeros(1,t_end);
flux_eis = zeros(1,t_end);
flux_rs  = zeros(1,t_end);
flux_rn  = zeros(1,t_end);
flux_esu = zeros(1,t_end);
flux_re  = zeros(1,t_end);
flux_qin = zeros(1,t_end);
flux_esb = zeros(1,t_end);
flux_qb  = zeros(1,t_end);
flux_qh  = zeros(1,t_end);
flux_qc  = zeros(1,t_end);
flux_ec  = zeros(1,t_end);
flux_qt  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Canopy storage
% S2. Stem/trunk storage
% S3. Upper soil layer
% S4. Lower soil layer
% S5. Hillslope routing
% S6. Channel routing

% RC(c,P(t)): rainfall on channel
RC = split_1;

% RG(1-c,P(t)): rainfall on ground
RG = split_1;

% EIC(S1,(1-c)*Ep(t),delta_t): evaporation from canopy
EIC = evap_1;

% QIE(RG(1-c,P(t)),S1,i1max): overflow from canopy
QIE = interception_1;

% RT(1-a,QIE(RG(1-c,P(t)),S1,i1max)): throughfall
RT = split_1;

% QIS(a,QIE(RG(1-c,P(t)),S1,i1max)): stem/trunk interception
QIS = split_1;

% EIS(S2,(1-c)*Ep(t),delta_t): evaporation from stem/trunk storage
EIS = evap_1;

% RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max): overflow from stem/trunk
RS = interception_1;

% RE(d50,d16,S3,RT(1-a,QIE(RG(1-c,P(t)),S1,i1max))+RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max)): 
% effective rain from surface
RE = saturation_13;

% ESU(S3,(1-c)*Ep(t),delta_t): evaporation from upper soil moisture
ESU = evap_1;

% QIN(kin,S3): infiltration 
QIN = recharge_3;

% ESB(1,S4,sbc,max(0,(1-c)*Ep(t)-ESU(S3,(1-c)*Ep(t)),delta_t): evaporation from deeper soil
ESB = evap_3;

% QB(kb,pb,S4,delta_t): baseflow
QB = baseflow_7;

% QH(kh,5/3,S5,delta_t): hillslope flow
QH = interflow_3;

% QC(kc,5/3,S6,delta_t): channel flow
QC = interflow_3;

% QT(QB(kb,pb,S4)+QH(kh,5/3,S5)+QC(kc,5/3,S6),c*Ep(t)): evaporation from channel
QT = effective_1;

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
                                               0,0,1,1,0,0;
                                               1,1,1,0,1,0;
                                               0,0,0,0,0,1]);               % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0,0,0;
                                                  1,1,0,0,0,0;
                                                  1,1,1,0,0,0;
                                                  0,0,1,1,0,0;
                                                  1,1,1,0,1,0;
                                                  0,0,0,0,0,1],...
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
                (RG(1-c,P(t)) - ...
                 EIC(S1,(1-c)*Ep(t),delta_t) - ...
                 QIE(RG(1-c,P(t)),S1,i1max));                               % Store 1 function with current flux values
    
    tmpf_S2 = @(S1,S2,S3,S4,S5,S6) ...
                (QIS(a,QIE(RG(1-c,P(t)),S1,i1max)) - ...
                 EIS(S2,(1-c)*Ep(t),delta_t) - ...
                 RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max));           % Store 2 function
            
    tmpf_S3 = @(S1,S2,S3,S4,S5,S6) ...
                (RT(1-a,QIE(RG(1-c,P(t)),S1,i1max)) + ...
                 RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max) - ...
                 ESU(S3,(1-c)*Ep(t),delta_t) - ...
                 RE(d50,d16,S3,RT(1-a,QIE(RG(1-c,P(t)),S1,i1max))+RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max)) - ...
                 QIN(kin,S3));                                              % Store 3 function
            
    tmpf_S4 = @(S1,S2,S3,S4,S5,S6) ...
                (QIN(kin,S3) - ...
                 ESB(1,S4,sbc,max(0,(1-c)*Ep(t)-ESU(S3,(1-c)*Ep(t),delta_t)),delta_t) - ...
                 QB(kb,pb,S4,delta_t));                                     % Store 4 function
            
    tmpf_S5 = @(S1,S2,S3,S4,S5,S6) ...
                (RE(d50,d16,S3,RT(1-a,QIE(RG(1-c,P(t)),S1,i1max))+RS(QIS(a,QIE(RG(1-c,P(t)),S1,i1max)),S2,i2max)) - ...
                 QH(kh,5/3,S5,delta_t));                                    % Store 5 function
            
    tmpf_S6 = @(S1,S2,S3,S4,S5,S6) ...
                (RC(c,P(t)) - ...
                 QC(kc,5/3,S6,delta_t));                                    % Store 6 function
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old,S6old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,tmpf_S5,tmpf_S6);        % this returns a new anonymous function that we solve in the next step

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
                                        eq_sys(4),eq_sys(5),eq_sys(6)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,...
                                         S4old,S5old,S6old], ...            % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_rc(t)  = RC(c,P(t));
    flux_rg(t)  = RG(1-c,P(t));
    flux_eic(t) = EIC(tmp_sFlux(1),(1-c)*Ep(t),delta_t);
    flux_qie(t) = QIE(flux_rg(t),tmp_sFlux(1),i1max);
    flux_qis(t) = QIS(a,flux_qie(t));
    flux_rt(t)  = RT(1-a,flux_qie(t));
    flux_eis(t) = EIS(tmp_sFlux(2),(1-c)*Ep(t),delta_t);
    flux_rs(t)  = RS(flux_qis(t),tmp_sFlux(2),i2max);
    flux_rn(t)  = flux_rt(t) + flux_rs(t);
    flux_esu(t) = ESU(tmp_sFlux(3),(1-c)*Ep(t),delta_t);
    flux_re(t)  = RE(d50,d16,tmp_sFlux(3),flux_rn(t));
    flux_qin(t) = QIN(kin,tmp_sFlux(3));
    flux_esb(t) = ESB(1,tmp_sFlux(4),sbc,max(0,(1-c)*Ep(t)-flux_esu(t)),delta_t);
    flux_qb(t)  = QB(kb,pb,tmp_sFlux(4),delta_t);
    flux_qh(t)  = QH(kh,5/3,tmp_sFlux(5),delta_t);
    flux_qc(t)  = QC(kc,5/3,tmp_sFlux(6),delta_t);
    flux_qt(t)  = QT(flux_qb(t)+flux_qh(t)+flux_qc(t),c*Ep(t));
    flux_ec(t)  = (flux_qb(t)+flux_qh(t)+flux_qc(t)) - flux_qt(t);
    
    % Update the stores
    store_S1(t) = S1old + (flux_rg(t)  - flux_eic(t) - flux_qie(t)) * delta_t;
    store_S2(t) = S2old + (flux_qis(t) - flux_eis(t) - flux_rs(t)) * delta_t;
    store_S3(t) = S3old + (flux_rn(t)  - flux_re(t)  - flux_esu(t) - flux_qin(t)) * delta_t;
    store_S4(t) = S4old + (flux_qin(t) - flux_esb(t) - flux_qb(t)) * delta_t;
    store_S5(t) = S5old + (flux_re(t)  - flux_qh(t)) * delta_t;
    store_S6(t) = S6old + (flux_rc(t)  - flux_qc(t)) * delta_t;

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eic + flux_eis + flux_esu + flux_esb +flux_ec) * delta_t;
    fluxOutput.Q      = (flux_qt) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.rc   = flux_rs * delta_t;
    fluxInternal.rg   = flux_rg * delta_t;
    fluxInternal.eic  = flux_eic * delta_t;
    fluxInternal.qie  = flux_qie * delta_t;
    fluxInternal.qis  = flux_qis * delta_t;
    fluxInternal.eis  = flux_eis * delta_t;
    fluxInternal.rs   = flux_rs * delta_t;
    fluxInternal.rn   = flux_rn * delta_t;
    fluxInternal.esu  = flux_esu * delta_t;
    fluxInternal.qin  = flux_qin * delta_t;
    fluxInternal.esb  = flux_esb * delta_t;
    fluxInternal.qb   = flux_qb * delta_t;
    fluxInternal.qh   = flux_qh * delta_t;
    fluxInternal.qc   = flux_qc * delta_t;
    fluxInternal.qt   = flux_qt * delta_t;
    fluxInternal.ec   = flux_ec * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P,...              % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


 




