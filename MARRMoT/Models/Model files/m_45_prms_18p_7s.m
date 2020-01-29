function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_45_prms_18p_7s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: PRMS
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Leavesley, G. H., R. Lichty, B. Troutman, and L. Saindon (1983), 
% Precipitation-Runo Modeling System: User's Manual. U.S. Geological 
% Survey, Water-Resources Investigations Report 83-4238, 207
%
% Markstrom, S. L., S. Regan, L. E. Hay, R. J. Viger, R. M. T. Webb, R. A. 
% Payn, and J. H. LaFontaine (2015), PRMS-IV, the Precipitation-Runoff 
% Modeling System, Version 4. In U.S. Geological Survey Techniques and
% Methods, book 6, chap. B7, 158. doi: http://dx.doi.org/10.3133/tm6B7

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
tt      = theta(1);     % Temperature threshold for snowfall and melt [oC]
ddf     = theta(2);     % Degree-day factor for snowmelt [mm/oC/d]
alpha   = theta(3);     % Fraction of rainfall on soil moisture going to interception [-] 
beta    = theta(4);     % Fraction of catchment where rain goes to soil moisture [-]
stor    = theta(5);     % Maximum interception capcity [mm]
retip   = theta(6);     % Maximum impervious area storage [mm]
fscn    = theta(7);     % Fraction of SCX where SCN is located [-]
scx     = theta(8);     % Maximum contributing fraction area to saturation excess flow [-]
scn     = fscn*scx;     % Minimum contributing fraction area to saturation excess flow [-]
flz     = theta(9);     % Fraction of total soil moisture that is the lower zone [-]
stot    = theta(10);    % Total soil moisture storage [mm]: REMX+SMAX
remx    = (1-flz)*stot; % Maximum upper soil moisture storage [mm]
smax    = flz*stot;     % Maximum lower soil moisture storage [mm] 
cgw     = theta(11);    % Constant drainage to deep groundwater [mm/d]
resmax  = theta(12);    % Maximum flow routing reservoir storage (used for scaling only, there is no overflow) [mm]
k1      = theta(13);    % Groundwater drainage coefficient [d-1]
k2      = theta(14);    % Groundwater drainage non-linearity [-]
k3      = theta(15);    % Interflow coefficient 1 [d-1]
k4      = theta(16);    % Interflow coefficient 2 [mm-1 d-1]
k5      = theta(17);    % Baseflow coefficient [d-1]
k6      = theta(18);    % Groundwater sink coefficient [d-1]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial snow pack
S20     = storeInitial(2);       % Initial interception storage
S30     = storeInitial(3);       % Initial impervious area storage
S40     = storeInitial(4);       % Initial upper soil moisture
S50     = storeInitial(5);       % Initial lower soil moisture
S60     = storeInitial(6);       % Initial runoff routing reservoir
S70     = storeInitial(7);       % Initial groundwater

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0,0];     % lower bounds of stores
store_upp = [];                  % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);
store_S7 = zeros(1,t_end);

flux_ps  = zeros(1,t_end);
flux_pr  = zeros(1,t_end);
flux_pim = zeros(1,t_end);
flux_psm = zeros(1,t_end);
flux_pby = zeros(1,t_end);
flux_pin = zeros(1,t_end);
flux_ptf = zeros(1,t_end);
flux_m   = zeros(1,t_end);
flux_mim = zeros(1,t_end);
flux_msm = zeros(1,t_end);
flux_sas = zeros(1,t_end);
flux_sro = zeros(1,t_end);
flux_inf = zeros(1,t_end);
flux_pc  = zeros(1,t_end);
flux_excs= zeros(1,t_end);
flux_qres= zeros(1,t_end);
flux_sep = zeros(1,t_end);
flux_gad = zeros(1,t_end);
flux_ras = zeros(1,t_end);
flux_bas = zeros(1,t_end);
flux_snk = zeros(1,t_end);
flux_ein = zeros(1,t_end);
flux_eim = zeros(1,t_end);
flux_ea  = zeros(1,t_end);
flux_et  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. snow pack
% S2. interception storage
% S3. impervious area storage
% S4. upper soil moisture
% S5. lower soil moisture
% S6. runoff routing reservoir
% S7. groundwater

% PS(P(t),T(t),tt): precip as snow
PS = snowfall_1;

% M(ddf,tt,T(t),S1,delta_t): snowmelt
M = melt_1;

% PR(P(t),T(t),tt): precip as rain
PR = rainfall_1;

% PIM(1-beta,PR(P(t),T(t),tt)): fraction rain on impervious area
PIM = split_1;

% PSM(beta,PR(P(t),T(t),tt)): fraction rain on infiltration area
PSM = split_1;

% PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))): fraction rain on infiltration
% area that bypasses interception
PBY = split_1;

% PIN(alpha,PSM(beta,PR(P(t),T(t),tt))): fraction rain on infiltration area
% going to interception
PIN = split_1;

% PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor): overflow from
% interception
PTF = interception_1;

% EIN(S2,beta*Ep(t),delta_t): evaporation from interception
EIN = evap_1;

% MIM(1-beta,M(ddf,tt,T(t),S1,delta_t)): fraction snowmelt on impervious area
MIM = split_1;

% EIM(S3,(1-beta)*Ep(t),delta_t): evaporation from impervious area
EIM = evap_1;

% SAS(PIM(1-beta,PR(P(t),T(t),tt))+MIM(1-beta,M(ddf,tt,T(t),S1,delta_t)),S3,retip): 
% surface runoff from impervious area
SAS = saturation_1;

% MSM(beta,M(ddf,tt,T(t),S1,delta_t)): fraction snowmelt on infiltration area
MSM = split_1;

% SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) ):
% surface runoff from variable contributing area
SRO = saturation_8;

% INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )): 
% infiltration after surface runoff
INF = effective_1;

% EA(S4,remx,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t),delta_t): evaporation from upper soil moisture
EA = evap_7;

% PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx): 
% percolation to lower soil moisture
PC = saturation_1;

% ET(Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t)-EA(S4,remx,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t)),S5,smax,S4,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t),delta_t): 
% transpiration from lower soil moisture
ET = evap_15;

% EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax): 
% excess flow to deeper groundwater
EXCS = saturation_1;

% SEP(cgw,EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax)): 
% minimum of constant drainage and flux availability
SEP = recharge_7;

% QRES(EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))))),S4,remx),S5,smax),SEP(cgw,EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax))):
% flow to runoff reservoir
QRES = effective_1;

% GAD(k2,S6,resmax,k1): groundwater recharge
GAD = recharge_2;

% RAS(k3,k4,S6): interflow
RAS = interflow_4;

% BAS(k5,S7): baseflow
BAS = baseflow_1;

% SNK(k6,S7): sink flow
SNK = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0,0;
                                               0,1,0,0,0,0,0;
                                               1,0,1,0,0,0,0;
                                               1,1,1,1,0,0,0;
                                               1,1,1,1,1,0,0;
                                               1,1,1,1,1,1,0;
                                               1,1,1,1,1,1,1]);             % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [  1,0,0,0,0,0,0;
                                                    0,1,0,0,0,0,0;
                                                    1,0,1,0,0,0,0;
                                                    1,1,1,1,0,0,0;
                                                    1,1,1,1,1,0,0;
                                                    1,1,1,1,1,1,0;
                                                    1,1,1,1,1,1,1],...
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
    if t == 1; S7old = S70; else; S7old = store_S7(t-1); end                % store 7 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 2 function
                (PS(P(t),T(t),tt) - ...
                 M(ddf,tt,T(t),S1,delta_t));                                

    tmpf_S2 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 2 function
                (PIN(alpha,PSM(beta,PR(P(t),T(t),tt))) - ...
                 EIN(S2,beta*Ep(t),delta_t) - ...
                 PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor));       

    tmpf_S3 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 3 function
                (PIM(1-beta,PR(P(t),T(t),tt)) + ...
                 MIM(1-beta,M(ddf,tt,T(t),S1,delta_t)) - ...
                 EIM(S3,(1-beta)*Ep(t),delta_t) - ...
                 SAS(PIM(1-beta,PR(P(t),T(t),tt))+MIM(1-beta,M(ddf,tt,T(t),S1,delta_t)),S3,retip));

    tmpf_S4 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 4 function
                (INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )) - ...
                 EA(S4,remx,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t),delta_t) - ...
                 PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx));

    tmpf_S5 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 5 function
                (PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx) - ...
                 ET(Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t)-EA(S4,remx,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t),delta_t),S5,smax,S4,Ep(t)-EIN(S2,beta*Ep(t),delta_t)-EIM(S3,(1-beta)*Ep(t),delta_t),delta_t) - ...
                 EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax));
    
    tmpf_S6 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 6 function
                (QRES(EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))))),S4,remx),S5,smax),SEP(cgw,EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax))) - ...
                 GAD(k2,S6,resmax,k1) - ...
                 RAS(k3,k4,S6));

    tmpf_S7 = @(S1,S2,S3,S4,S5,S6,S7) ...                                   % Store 7 function
                (SEP(cgw,EXCS(PC(INF(MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))),SRO(scn,scx,S4,remx,MSM(beta,M(ddf,tt,T(t),S1,delta_t))+PTF(PIN(alpha,PSM(beta,PR(P(t),T(t),tt))),S2,stor)+PBY(1-alpha,PSM(beta,PR(P(t),T(t),tt))) )),S4,remx),S5,smax)) + ...
                 GAD(k2,S6,resmax,k1) - ...
                 BAS(k5,S7) - ...
                 SNK(k6,S7));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,S5old,S6old,S7old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,...
                        tmpf_S5,tmpf_S6,tmpf_S7);                           % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                            eq_sys(5),eq_sys(6),eq_sys(7)),...              % system of storage equations
                        [S1old,S2old,S3old,S4old,S5old,S6old,S7old],...                                   % storage values on previous time step
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
                                        eq_sys(4),eq_sys(5),eq_sys(6),...
                                        eq_sys(7)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old,...
                                         S5old,S6old,S7old], ...            % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ps(t)  = PS(P(t),T(t),tt);
    flux_pr(t)  = PR(P(t),T(t),tt);
    flux_pim(t) = PIM(1-beta,flux_pr(t));
    flux_psm(t) = PSM(beta,flux_pr(t));
    flux_pby(t) = PBY(1-alpha,flux_psm(t));
    flux_pin(t) = PIN(alpha,flux_psm(t));
    flux_ptf(t) = PTF(flux_pin(t),tmp_sFlux(2),stor);
    flux_m(t)   = M(ddf,tt,T(t),tmp_sFlux(1),delta_t);
    flux_mim(t) = MIM(1-beta,flux_m(t));
    flux_msm(t) = MSM(beta,flux_m(t));
    flux_sas(t) = SAS(flux_pim(t)+flux_mim(t),tmp_sFlux(3),retip);
    flux_sro(t) = SRO(scn,scx,tmp_sFlux(4),remx,flux_msm(t)+flux_ptf(t)+flux_pby(t));
    flux_inf(t) = INF(flux_msm(t)+flux_ptf(t)+flux_pby(t),flux_sro(t));
    flux_pc(t)  = PC(flux_inf(t),tmp_sFlux(4),remx);
    flux_excs(t)= EXCS(flux_pc(t),tmp_sFlux(5),smax);
    flux_sep(t) = SEP(cgw,flux_excs(t));
    flux_qres(t)= QRES(flux_excs(t),flux_sep(t));
    flux_gad(t) = GAD(k2,tmp_sFlux(6),resmax,k1);
    flux_ras(t) = RAS(k3,k4,tmp_sFlux(6));
    flux_bas(t) = BAS(k5,tmp_sFlux(7));
    flux_snk(t) = SNK(k6,tmp_sFlux(7));
    flux_ein(t) = EIN(tmp_sFlux(2),beta*Ep(t),delta_t);
    flux_eim(t) = EIM(tmp_sFlux(3),(1-beta)*Ep(t),delta_t);
    flux_ea(t)  = EA(tmp_sFlux(4),remx,Ep(t)-flux_ein(t)-flux_eim(t),delta_t);
    flux_et(t)  = ET(Ep(t)-flux_ein(t)-flux_eim(t)-flux_ea(t),tmp_sFlux(5),smax,tmp_sFlux(4),Ep(t)-flux_ein(t)-flux_eim(t),delta_t);
        
    % Update the stores
    store_S1(t) = S1old + (flux_ps(t)  - flux_m(t)) * delta_t;
    store_S2(t) = S2old + (flux_pin(t) - flux_ein(t) - flux_ptf(t)) * delta_t;    
    store_S3(t) = S3old + (flux_pim(t) + flux_mim(t) - flux_eim(t) - flux_sas(t)) * delta_t;
    store_S4(t) = S4old + (flux_inf(t) - flux_ea(t)  - flux_pc(t)) * delta_t;
    store_S5(t) = S5old + (flux_pc(t)  - flux_et(t)  - flux_excs(t)) * delta_t;
    store_S6(t) = S6old + (flux_qres(t)- flux_gad(t) - flux_ras(t)) * delta_t;
    store_S7(t) = S7old + (flux_sep(t) + flux_gad(t) - flux_bas(t) - flux_snk(t)) * delta_t;
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ein + flux_eim + flux_ea + flux_et) * delta_t;
    fluxOutput.Q      = (flux_sas + flux_sro + flux_ras + flux_bas) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ps   = flux_ps * delta_t;
    fluxInternal.pr   = flux_pr * delta_t;
    fluxInternal.pim  = flux_pim * delta_t;
    fluxInternal.psm  = flux_psm * delta_t;
    fluxInternal.pby  = flux_pby * delta_t;
    fluxInternal.pin  = flux_pin * delta_t;
    fluxInternal.ptf  = flux_ptf * delta_t;
    fluxInternal.m    = flux_m * delta_t;
    fluxInternal.mim  = flux_mim * delta_t;
    fluxInternal.msm  = flux_msm * delta_t;
    fluxInternal.sas  = flux_sas * delta_t;
    fluxInternal.sro  = flux_sro * delta_t;
    fluxInternal.inf  = flux_inf * delta_t;
    fluxInternal.pc   = flux_pc * delta_t;
    fluxInternal.excs = flux_excs * delta_t;
    fluxInternal.sep  = flux_sep * delta_t;
    fluxInternal.qres = flux_qres * delta_t;
    fluxInternal.gad  = flux_gad * delta_t;
    fluxInternal.ras  = flux_ras * delta_t;
    fluxInternal.bas  = flux_bas * delta_t;
    fluxInternal.snk  = flux_snk * delta_t;
    fluxInternal.ein  = flux_ein * delta_t;
    fluxInternal.eim  = flux_eim * delta_t;
    fluxInternal.ea   = flux_ea * delta_t;
    fluxInternal.et   = flux_et * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;
    storeInternal.S7  = store_S7;

% Check water balance
if nargout == 4

    % Manual water balance due to sink flow
    waterBalance =  sum(P).*delta_t - ...       % Incoming precipitation
                    sum(fluxOutput.Q) - ...     % Outgoing flow
                    sum(fluxOutput.Ea) - ...    % Outgoing evaporation
                    sum(flux_snk) - ...         % Outgoing sink flow
                    (store_S1(end)-S10) - ...   % Storage change
                    (store_S2(end)-S20) - ...
                    (store_S3(end)-S30) - ...
                    (store_S4(end)-S40) - ...
                    (store_S5(end)-S50) - ...
                    (store_S6(end)-S60) - ...
                    (store_S7(end)-S70);

    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
    disp(['Delta S4 = ',num2str((store_S4(end)-S40)),'mm.'])
    disp(['Delta S5 = ',num2str((store_S5(end)-S50)),'mm.'])
    disp(['Delta S6 = ',num2str((store_S6(end)-S60)),'mm.'])
    disp(['Delta S7 = ',num2str((store_S7(end)-S70)),'mm.'])
    disp(['Sink flow= ',num2str(sum(flux_snk)),'mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(sink)) - delta S = ',num2str(waterBalance)]);     

end




