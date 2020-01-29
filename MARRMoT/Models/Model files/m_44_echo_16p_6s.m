function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_44_echo_16p_6s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: ECHO
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
%
% Schaefli, B., Nicotina, L., Imfeld, C., Da Ronco, P., Bertuzzo, E., & 
% Rinaldo, A. (2014). SEHR-ECHO v1.0: A spatially explicit hydrologic 
% response model for ecohydrologic applications. Geoscientific Model 
% Development, 7(6), 2733–2746. http://doi.org/10.5194/gmd-7-2733-2014


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
rho  = theta(1);     % Maximum interception storage [mm]
ts   = theta(2);     % Threshold temperature for snowfall [oC]
tm   = theta(3);     % Threshold temperature for snowmelt [oC]
as   = theta(4);     % Degree-day factor [mm/oC/d]
af   = theta(5);     % Refreezing reduction factor [-] 
gmax = theta(6);     % Maximum melt due to ground-heat flux [mm/d]
the  = theta(7);     % Water-holding capacity of snow [-]
phi  = theta(8);     % Maximum infiltration rate [mm/d]
smax = theta(9);     % Maximum soil moisture storage [mm]
fsm  = theta(10);    % Plant stress point as a fraction of Smax [-]
fsw  = theta(11);    % Wilting point as fraction of sm [-]
ksat = theta(12);    % Runoff rate from soil moisture [d-1]
c    = theta(13);    % Runoff non-linearity from soil moisture [-]
lmax = theta(14);    % Groundwater flux [mm/d]
kf   = theta(15);    % Runoff coefficient [d-1]
ks   = theta(16);    % Runoff coefficient [d-1]

% Auxiliary parameters
sm   = fsm*smax;     % Plant stress point [mm]
sw   = fsw*sm;       % Wilting point [mm]

%%INITIALISE MODEL STORES
S10  = storeInitial(1);       % Initial interception storage
S20  = storeInitial(2);       % Initial snow storage
S30  = storeInitial(3);       % Initial snow water storage
S40  = storeInitial(4);       % Initial soil moisture storage
S50  = storeInitial(5);       % Initial fast routing storage
S60  = storeInitial(6);       % Initial slow routing storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0];    % lower bounds of stores
store_upp = [];               % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);

flux_ei  = zeros(1,t_end);
flux_pn  = zeros(1,t_end);
flux_ps  = zeros(1,t_end);
flux_pr  = zeros(1,t_end);
flux_ms  = zeros(1,t_end);
flux_fs  = zeros(1,t_end);
flux_gs  = zeros(1,t_end);
flux_mw  = zeros(1,t_end);
flux_ew  = zeros(1,t_end);
flux_eq  = zeros(1,t_end);
flux_rh  = zeros(1,t_end);
flux_eps = zeros(1,t_end);
flux_et  = zeros(1,t_end);
flux_fi  = zeros(1,t_end);
flux_rd  = zeros(1,t_end);
flux_l   = zeros(1,t_end);
flux_lf  = zeros(1,t_end);
flux_ls  = zeros(1,t_end);
flux_rf  = zeros(1,t_end);
flux_rs  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception
% S2. Snow
% S3. Water in snow
% S4. Soil moisture
% S5. Fast routing
% S6. Slow routing

% EI(S1,Ep(t),delta_t): evaporation from interception
EI = evap_1;

% PN(P(t),S1,rho): overflow from interception
PN = interception_1;

% PS(PN(P(t),S1,rho),T(t),ts): snowfall
PS = snowfall_1;

% PR(PN(P(t),S1,rho),T(t),ts): rainfall
PR = rainfall_1;

% MS(as,tm,T(t),S2,delta_t): snow melt
MS = melt_1;

% FS(af,as,tm,T(t),S3,delta_t): refreezing of liquid water in snow
FS = refreeze_1;

% GS(gmax,S2,delta_t): melt due to ground-heat flux
GS = melt_2;

% MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2,delta_t),S3,the*S2): overflow 
% from snowpack storage
MW = saturation_1;

% EW(S3old,the*S2,delta_t): outflow of excess storage if store size changes 
EW = excess_1;

% FI(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2),S3,the*S2)+GS(gmax,S2)+EW(S3old,the*S2),phi): 
% infiltration
FI = infiltration_4;

% RH(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2),S3,the*S2)+GS(gmax,S2),FI(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2),S3,the*S2)+GS(gmax,S2),phi)): 
% infiltration excess runoff
RH = effective_1;

% RD(FI(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2),S3,the*S2)+GS(gmax,S2),phi),S4,smax):
% saturation excess
RD = saturation_1;

% EPS(Ep(t),EI(S1,Ep(t))): unfulfilled EP demand
EPS = effective_1;

% ET(sw,sm,S4,EPS(Ep(t),EI(S1,Ep(t),delta_t)),delta_t): evapotranspiration
ET = evap_22;

% L(ksat,c,S4,delta_t): leakage from soil moisture
L = recharge_6;

% LS(lmax,L(ksat,c,S4)): leakage to slow routing
LS = recharge_7;

% LF(L(ksat,c,S4),LS(lmax,L(ksat,c,S4))): leakage to fast routing
LF = effective_1;

% RF(kf,S5): fast flow
RF = baseflow_1;

% RS(ks,S6): slow flow
RS = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,0,0,0;
                                               1,1,1,1,0,0;
                                               0,0,0,1,1,0;
                                               0,0,0,1,0,1]);               % Specify the Jacobian pattern                                               

lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0,0,0;
                                                  1,1,1,0,0,0;
                                                  1,1,1,0,0,0;
                                                  1,1,1,1,0,0;
                                                  0,0,0,1,1,0;
                                                  0,0,0,1,0,1],...
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
               (P(t) - ...
                EI(S1,Ep(t),delta_t) - ...
                PN(P(t),S1,rho));
            
	tmpf_S2 = @(S1,S2,S3,S4,S5,S6) ...
               (PS(PN(P(t),S1,rho),T(t),ts) +...
                FS(af,as,tm,T(t),S3,delta_t) - ...
                MS(as,tm,T(t),S2,delta_t) - ...
                GS(gmax,S2,delta_t));

    tmpf_S3 = @(S1,S2,S3,S4,S5,S6) ...
               (PR(PN(P(t),S1,rho),T(t),ts) + ...
                MS(as,tm,T(t),S2,delta_t) - ...
                FS(af,as,tm,T(t),S3,delta_t) - ...
                MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2,delta_t),S3,the*S2) - ...
                EW(S3old,the*S2,delta_t));

    tmpf_S4 = @(S1,S2,S3,S4,S5,S6) ...
               (FI(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2,delta_t),S3,the*S2)+GS(gmax,S2,delta_t)+EW(S3old,the*S2,delta_t),phi) - ...
                ET(sw,sm,S4,EPS(Ep(t),EI(S1,Ep(t),delta_t)),delta_t) - ...
                RD(FI(MW(PR(PN(P(t),S1,rho),T(t),ts)+MS(as,tm,T(t),S2,delta_t),S3,the*S2)+GS(gmax,S2,delta_t),phi),S4,smax) - ...
                L(ksat,c,S4,delta_t));

    tmpf_S5 = @(S1,S2,S3,S4,S5,S6) ...
               (LF(L(ksat,c,S4,delta_t),LS(lmax,L(ksat,c,S4,delta_t))) - ...
                RF(kf,S5));

    tmpf_S6 = @(S1,S2,S3,S4,S5,S6) ...
               (LS(lmax,L(ksat,c,S4,delta_t)) - ...
                RS(ks,S6));
    
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
                                        eq_sys(1),eq_sys(2),...
                                        eq_sys(3),eq_sys(4),...
                                        eq_sys(5),eq_sys(6)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,...
                                        S4old,S5old,S6old], ...             % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ei(t) = EI(tmp_sFlux(1),Ep(t),delta_t);
    flux_pn(t) = PN(P(t),tmp_sFlux(1),rho);
    flux_ps(t) = PS(flux_pn(t),T(t),ts);
    flux_pr(t) = PR(flux_pn(t),T(t),ts);
    flux_ms(t) = MS(as,tm,T(t),tmp_sFlux(2),delta_t);
    flux_fs(t) = FS(af,as,tm,T(t),tmp_sFlux(3),delta_t);
    flux_gs(t) = GS(gmax,tmp_sFlux(2),delta_t);
    flux_mw(t) = MW(flux_pr(t)+flux_ms(t),tmp_sFlux(3),the*tmp_sFlux(2));
    flux_ew(t) = EW(S3old,the*tmp_sFlux(2),delta_t);
    flux_eq(t) = flux_mw(t) + flux_gs(t) + flux_ew(t);
    flux_fi(t) = FI(flux_eq(t),phi);
    flux_rh(t) = RH(flux_eq(t),flux_fi(t));
    flux_eps(t)= EPS(Ep(t),flux_ei(t));
    flux_et(t) = ET(sw,sm,tmp_sFlux(4),flux_eps(t),delta_t);
    flux_rd(t) = RD(flux_fi(t),tmp_sFlux(4),smax);
    flux_l(t)  = L(ksat,c,tmp_sFlux(4),delta_t);
    flux_ls(t) = LS(lmax,flux_l(t));
    flux_lf(t) = LF(flux_l(t),flux_ls(t));
    flux_rf(t) = RF(kf,tmp_sFlux(5));
    flux_rs(t) = RS(ks,tmp_sFlux(6));
    
    % Update the stores
    store_S1(t) = S1old + (P(t)       - flux_ei(t) - flux_pn(t)) * delta_t;
    store_S2(t) = S2old + (flux_ps(t) + flux_fs(t) - flux_ms(t) - flux_gs(t)) * delta_t;
    store_S3(t) = S3old + (flux_pr(t) + flux_ms(t) - flux_fs(t) - flux_mw(t) - flux_ew(t)) * delta_t;
    store_S4(t) = S4old + (flux_fi(t) - flux_et(t) - flux_rd(t) - flux_l(t)) * delta_t;
    store_S5(t) = S5old + (flux_lf(t) - flux_rf(t)) * delta_t;
    store_S6(t) = S6old + (flux_ls(t) - flux_rs(t)) * delta_t;

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ei + flux_et) * delta_t;
    fluxOutput.Q      = (flux_rh + flux_rd + flux_rf + flux_rs) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ei  = flux_ei * delta_t;
    fluxInternal.pn  = flux_pn * delta_t;
    fluxInternal.ps  = flux_ps * delta_t;
    fluxInternal.pr  = flux_pr * delta_t;
    fluxInternal.ms  = flux_ms * delta_t;
    fluxInternal.fs  = flux_fs * delta_t;
    fluxInternal.gs  = flux_gs * delta_t;
    fluxInternal.mw  = flux_mw * delta_t;
    fluxInternal.ew  = flux_ew * delta_t;
    fluxInternal.eq  = flux_eq * delta_t;
    fluxInternal.fi  = flux_fi * delta_t;
    fluxInternal.rh  = flux_rh * delta_t;
    fluxInternal.eps = flux_eps * delta_t;
    fluxInternal.et  = flux_et * delta_t;
    fluxInternal.rd  = flux_rd * delta_t;
    fluxInternal.l   = flux_l * delta_t;
    fluxInternal.ls  = flux_ls * delta_t;
    fluxInternal.lf  = flux_lf * delta_t;
    fluxInternal.rf  = flux_rf * delta_t;
    fluxInternal.rs  = flux_rs * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


 




