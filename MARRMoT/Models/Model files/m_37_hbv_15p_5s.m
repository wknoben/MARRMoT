function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
                m_37_hbv_15p_5s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: HBV-96
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Lindström, G., Johansson, B., Persson, M., Gardelin, M., & Bergström, S. 
% (1997). Development and test of the distributed HBV-96 hydrological model. 
% Journal of Hydrology, 201, 272–288. 
% https://doi.org/https://doi.org/10.1016/S0022-1694(97)00041-3


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
tt      = theta(1);     % TT, middle of snow-rain interval [oC]
tti     = theta(2);     % TTI, interval length of rain-snow spectrum [oC]
ttm     = theta(3);     % TTM, threshold temperature for snowmelt [oC]
cfr     = theta(4);     % CFR, coefficient of refreezing of melted snow [-]
cfmax   = theta(5);     % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
whc     = theta(6);     % WHC, maximum water holding content of snow pack [-]
cflux   = theta(7);   	% CFLUX, maximum rate of capillary rise [mm/d]
fc      = theta(8);     % FC, maximum soil moisture storage [mm]
lp      = theta(9);     % LP, wilting point as fraction of FC [-]
beta    = theta(10);    % BETA, non-linearity coefficient of upper zone recharge [-]
k0      = theta(11);    % K0, runoff coefficient from upper zone [d-1], 
alpha   = theta(12);    % ALPHA, non-linearity coefficient of runoff from upper zone [-]
perc    = theta(13);    % PERC, maximum rate of percolation to lower zone [mm/d]
k1      = theta(14);    % K1, runoff coefficient from lower zone [d-1]
maxbas  = theta(15);    % MAXBAS, flow routing delay [d]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial snow pack storage
S20     = storeInitial(2);       % Initial water content of snow pack
S30     = storeInitial(3);       % Initial soil moisture storage
S40     = storeInitial(4);       % Initial upper zone storage
S50     = storeInitial(5);       % Initial lower zone storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0];           % lower bounds of stores
store_upp = [];                    % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);         % Snow pack
store_S2 = zeros(1,t_end);         % Water content of snow
store_S3 = zeros(1,t_end);         % Soil moisture
store_S4 = zeros(1,t_end);         % Upper zone
store_S5 = zeros(1,t_end);         % Lower zone

flux_sf     = zeros(1,t_end);         
flux_refr   = zeros(1,t_end);
flux_melt   = zeros(1,t_end);
flux_rf     = zeros(1,t_end);
flux_in     = zeros(1,t_end);
flux_se     = zeros(1,t_end);
flux_cf     = zeros(1,t_end);
flux_ea     = zeros(1,t_end);
flux_r      = zeros(1,t_end);
flux_q0     = zeros(1,t_end);
flux_perc   = zeros(1,t_end);
flux_q1     = zeros(1,t_end);
flux_qt     = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
% [Optional]
[~,uh_full] = uh_4_full(1,maxbas,delta_t);  % the Unit Hydrograph function spreads flow over maxbas days

%%INITIALISE ROUTING VECTORS
tmp_Qt_old = zeros(1,length(uh_full));      % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Snow pack
% S2. Water content of snow
% S3. Soil moisture
% S4. Upper zone
% S5. Lower zone

% SF(P,T,tt,tti): snowfall
SF = snowfall_2;

% REFR(cfr,cfmax,ttm,T,S2,delta_t): refreezing of melted snow
REFR = refreeze_1;

% MELT(cfmax,ttm,T,S1,delta_t): snowmelt
MELT = melt_1;

% RF(P,T,tt,tti): rainfall
RF = rainfall_2;

% IN(RF+MELT,S2,whc*S1): excess of liquid-water-in-snow storage
IN = infiltration_3;

% SE(S2old,whc*S1,delta_t): storage excess when store size changes
SE = excess_1;

% CF(cflux,S3,fc,S4,delta_t): capillary rise
CF = capillary_1;

% EA(lp,S3,fc,Ep,delta_t): evaporation from soil moisture
EA = evap_3;

% R(beta,S3,fc,IN): recharge to upper zone
R = recharge_2;

% Q0(k0,S4,alpha,delta_t): interflow
Q0 = interflow_2;

% PERC(perc,S4,delta_t): percolation to lower zone
PERC = percolation_1;

% Q1(k1,S5): baseflow
Q1 = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...
                                'JacobPattern', [1,1,0,0,0;
                                                 1,1,0,0,0;
                                                 1,1,1,1,0;
                                                 1,1,1,1,0;
                                                 0,0,0,1,1]);               % Specify the Jacobian pattern                                             
lsqnonlin_options = optimoptions('lsqnonlin',...
                                 'Display','none',...
                                'JacobPattern', [1,1,0,0,0;
                                                 1,1,0,0,0;
                                                 1,1,1,1,0;
                                                 1,1,1,1,0;
                                                 0,0,0,1,1],...
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
    tmpf_S1 = @(S1,S2,S3,S4,S5) (SF(P(t),T(t),tt,tti) + ...
                                 REFR(cfr,cfmax,ttm,T(t),S2,delta_t) - ...
                                 MELT(cfmax,ttm,T(t),S1,delta_t));          % store 1 function with current flux values

    tmpf_S2 = @(S1,S2,S3,S4,S5) (RF(P(t),T(t),tt,tti) + ...
                                 MELT(cfmax,ttm,T(t),S1,delta_t) - ...
                                 REFR(cfr,cfmax,ttm,T(t),S2,delta_t) - ...
                                 IN(RF(P(t),T(t),tt,tti) + ...
                                    MELT(cfmax,ttm,T(t),S1,delta_t),S2,whc*S1) - ...
                                 SE(S2old,whc*S1,delta_t));                 % Store 2 function

    tmpf_S3 = @(S1,S2,S3,S4,S5) (IN(RF(P(t),T(t),tt,tti) + ...
                                    MELT(cfmax,ttm,T(t),S1,delta_t),S2,whc*S1) + ...
                                 SE(S2old,whc*S1,delta_t) + ...
                                 CF(cflux,S3,fc,S4,delta_t) - ...
                                 EA(lp,S3,fc,Ep(t),delta_t) - ...
                                 R(beta,S3,fc,...
                                   IN(RF(P(t),T(t),tt,tti) + ...
                                      MELT(cfmax,ttm,T(t),S1,delta_t),S2,whc*S1) + ...
                                      SE(S2old,whc*S1,delta_t)));           % Store 3 function

    tmpf_S4 = @(S1,S2,S3,S4,S5) (R(beta,S3,fc,...
                                   IN(RF(P(t),T(t),tt,tti) + ...
                                      MELT(cfmax,ttm,T(t),S1,delta_t),S2,whc*S1) + ...
                                      SE(S2old,whc*S1,delta_t)) - ...
                                 CF(cflux,S3,fc,S4,delta_t) - ...
                                 Q0(k0,S4,alpha,delta_t) - ...
                                 PERC(perc,S4,delta_t));                    % Store 4 function

    tmpf_S5 = @(S1,S2,S3,S4,S5) (PERC(perc,S4,delta_t) - ...
                                 Q1(k1,S5));                                % Store 5 function

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
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...
                                            eq_sys(1),eq_sys(2),eq_sys(3),...   
                                            eq_sys(4),eq_sys(5)), ...       % system of ODEs
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
    flux_sf(t)   = SF(P(t),T(t),tt,tti);
    flux_refr(t) = REFR(cfr,cfmax,ttm,T(t),tmp_sFlux(2),delta_t);
    flux_melt(t) = MELT(cfmax,ttm,T(t),tmp_sFlux(1),delta_t);
    flux_rf(t)   = RF(P(t),T(t),tt,tti);
    flux_in(t)   = IN(flux_rf(t)+flux_melt(t),tmp_sFlux(2),whc*tmp_sFlux(1));
    flux_se(t)   = SE(S2old,whc*tmp_sFlux(1),delta_t);
    flux_cf(t)   = CF(cflux,tmp_sFlux(3),fc,tmp_sFlux(4),delta_t);
    flux_ea(t)   = EA(lp,tmp_sFlux(3),fc,Ep(t),delta_t);
    flux_r(t)    = R(beta,tmp_sFlux(3),fc,flux_in(t)+flux_se(t));
    flux_q0(t)   = Q0(k0,tmp_sFlux(4),alpha,delta_t);
    flux_perc(t) = PERC(perc,tmp_sFlux(4),delta_t);
    flux_q1(t)   = Q1(k1,tmp_sFlux(5));
    
    % Update the stores
    store_S1(t) = S1old + (flux_sf(t)   + flux_refr(t) - flux_melt(t)) * delta_t;
    store_S2(t) = S2old + (flux_rf(t)   + flux_melt(t) - flux_refr(t) - flux_in(t) - flux_se(t)) * delta_t;
    store_S3(t) = S3old + (flux_in(t)   + flux_se(t)   + flux_cf(t)   - flux_ea(t) - flux_r(t)) * delta_t;
    store_S4(t) = S4old + (flux_r(t)    - flux_cf(t)   - flux_q0(t)   - flux_perc(t)) * delta_t;
    store_S5(t) = S5old + (flux_perc(t) - flux_q1(t)) * delta_t;
    
    % Total runoff Qt = Q0 + Q1. Apply a triangular routing scheme with
    % time base 'maxbas'
    tmp_Qt_cur    = (flux_q0(t) + flux_q1(t)).*uh_full;                     % find how the runoff of this time step will be distributed in time
    tmp_Qt_old    = tmp_Qt_old + tmp_Qt_cur;                                % update the 'still-to-flow-out' vector
    flux_qt(t)    = tmp_Qt_old(1);                                          % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_Qt_old = circshift(tmp_Qt_old,-1);                                  % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_Qt_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.sf   = flux_sf * delta_t;
    fluxInternal.refr = flux_refr * delta_t;
    fluxInternal.melt = flux_melt * delta_t;
    fluxInternal.rf   = flux_rf * delta_t;
    fluxInternal.in   = flux_in * delta_t;
    fluxInternal.se   = flux_se * delta_t;
    fluxInternal.cf   = flux_cf * delta_t;
    fluxInternal.ea   = flux_ea * delta_t;
    fluxInternal.r    = flux_r * delta_t;
    fluxInternal.q0   = flux_q0 * delta_t;
    fluxInternal.perc = flux_perc * delta_t;
    fluxInternal.q1   = flux_q1 * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    
% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...         % Incoming precipitation
                                     fluxOutput,...         % Fluxes Q and Ea leaving the model
                                     storeInternal,...      % Time series of storages ...
                                     storeInitial,...       % And initial store values to calculate delta S
                                     tmp_Qt_old.*delta_t);  % Whether the model uses a routing scheme that
                                                            % still contains water. Use '0' for no routing
end                                                         






