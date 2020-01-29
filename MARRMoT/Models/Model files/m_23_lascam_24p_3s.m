function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_23_lascam_24p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Large-scale catchment water and salt balance
% model element
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Sivapalan, M., Ruprecht, J. K., & Viney, N. R. (1996). Water and salt 
% balance modelling to predict the effects of land-use changes in forested 
% catchments. 1. Small catchment water balance model. Hydrological 
% Processes, 10(3).


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

% Parameters (all lower case)
% [name in documentation] = theta(order in which specified in parameter file)
af   = theta(1);     % Catchment-scale infiltration parameter [mm/d]
bf   = theta(2);     % Catchment-scale infiltration non-linearity parameter [-]
stot = theta(3);     % Total catchment storage [mm]
xa   = theta(4);     % Fraction of Stot that is Amax [-]
xf   = theta(5);     % Fraction of Stot-Amx that is depth Fmax [-]
na   = theta(6);     % Fraction of Amax that is Amin [-];
ac   = theta(7);     % Variable contributing area scaling [-]
bc   = theta(8);     % Variable contributing area non-linearity [-]
ass  = theta(9);     % Subsurface saturation area scaling [-]
bss  = theta(10);    % Subsurface saturation area non-linearity [-]
c    = theta(11);    % Maximum infiltration rate [mm/d]
ag   = theta(12);    % Interception base parameter [mm/d]
bg   = theta(13);    % Interception fraction parameter [-]
gf   = theta(14);    % F-store evaporation scaling [-]
df   = theta(15);    % F-store evaporation non-linearity [-]
td   = theta(16);    % Recharge time parameter [d-1]
ab   = theta(17);    % Groundwater flow scaling [-]
bb   = theta(18);    % Groundwater flow base rate [mm/d]
ga   = theta(19);    % A-store evaporation scaling [-]
da   = theta(20);    % A-store evaporation non-linearity [-]
aa   = theta(21);    % Subsurface storm flow rate [mm/d]
ba   = theta(22);    % Subsurface storm flow non-linearity [-]
gb   = theta(23);    % B-store evaporation scaling [-]
db   = theta(24);    % B-store evaporation non-linearity [-]

% Auxiliary parameters
amax = xa*stot;             % Maximum contributing area depth [mm]
fmax = xf*(stot-amax);      % Infiltration depth scaling [mm]
bmax = (1-xf)*(stot-amax);  % Groundwater depth scaling [mm]
amin = na*amax;             % Minimum contributing area depth [mm]

%%INITIALISE MODEL STORES
S10  = storeInitial(1);     % Initial infiltration storage
S20  = storeInitial(2);     % Initial soil moisture storage
S30  = storeInitial(3);     % Initial saturated zone storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];        % lower bounds of stores
store_upp = [];             % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1  = zeros(1,t_end);
store_S2  = zeros(1,t_end);
store_S3  = zeros(1,t_end);

flux_ei   = zeros(1,t_end);
flux_pg   = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qie  = zeros(1,t_end);
flux_pc   = zeros(1,t_end);
flux_qsse = zeros(1,t_end);
flux_qsie = zeros(1,t_end);
flux_fa   = zeros(1,t_end);
flux_ef   = zeros(1,t_end);
flux_rf   = zeros(1,t_end);
flux_ea1  = zeros(1,t_end); % this splits the original equation (Ea = phi_c*Ep + ga*(A/Amax)^da*Ep) into two parts
flux_ea2  = zeros(1,t_end); % first part: phi_c*Ep. second part: ga*(A/Amax)^da*Ep
flux_qa   = zeros(1,t_end);
flux_ra   = zeros(1,t_end);
flux_qb   = zeros(1,t_end);
flux_eb   = zeros(1,t_end);

tmp_phiss = zeros(1,t_end);
tmp_phic  = zeros(1,t_end);
tmp_fss   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Infiltration store (F)
% S2. Upper permeable zone (A)
% S3. Lower saturated zone (B)

% PG(bg,ag,P(t)): stylized interception
PG = interception_5;

% EI(P(t),PG(bg,ag,P(t))) evaporation assumed from interception
EI = effective_1;

% QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))): saturation excess flow
QSE = saturation_11;

% PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c): actual 
% infiltration into subsurface
PC = infiltration_4;

% QIE(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)):
% infiltration excess
QIE = effective_1;

% PHISS(ass,bss,S2,amin,amax): subsurface contriuting area
PHISS = area_1;

% PHIC(ac,bc,S2,amin,amax): surface contriuting area
PHIC = area_1;

% QSSE(PHISS(ass,bss,S2,amin,amax),PHIC(ac,bc,S2,amin,amax),PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)):
QSSE = saturation_12;

% FSS(af,bf,S3,bmax,S1,fmax): infiltration rate
FSS = infiltration_5;

% FA(max(0,PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)*min(1,(1-PHISS(ass,bss,S2,amin,amax))/(1-PHIC(ac,bc,S2,amin,amax)))),FSS(af,bf,S3,bmax,S1,fmax)): 
% infiltration after subsurface runoff
FA = infiltration_4;

% QSIE(PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c),FA(PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)*(1-PHISS(ass,bss,S2,amin,amax))/(1-PHIC(ac,bc,S2,amin,amax)),FSS(af,bf,S3,bmax,S1,fmax))+QSSE(PHISS(ass,bss,S2,amin,amax),PHIC(ac,bc,S2,amin,amax),PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c))):
QSIE = effective_1;

% EF(gf,df,S1,fmax,Ep(t),delta_t): evaporation from F
EF = evap_19;

% RF(td,S1): recharge to groundwater
RF = recharge_3;

% EA1(S2,PHIC(ac,bc,S2,amin,amax)*Ep(t),delta_t): evaporation from A (phi_c*Ep)
EA1 = evap_1;

% EA2(ga,da,S2,amax,Ep(t),delta_t): evaporation from A (ga*(A/Amax)^da*Ep)
EA2 = evap_19;

% QA(aa,ba,S2,amin,amax,1): subsurface storm flow (the '1' disables an
% unused parameter in the equation)
QA = saturation_11;

% RA(PHIC(ac,bc,S2,amin,amax)*FSS(af,bf,S3,bmax,S1,fmax),S2,delta_t): recharge to groundwater
RA = recharge_4;

% EB(gb,db,S3,bmax,Ep(t),delta_t): evaporation from B
EB = evap_19;

% QB(bb,ab,S3,bmax): outflow from B
QB = baseflow_8;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1,1;
                                               1,1,1;
                                               1,1,1]);                     % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,1,1;
                                                  1,1,1;
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
    tmpf_S1 = @(S1,S2,S3) ...
        (FA(max(0,PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)*min(1,(1-PHISS(ass,bss,S2,amin,amax))/(1-PHIC(ac,bc,S2,amin,amax)))),FSS(af,bf,S3,bmax,S1,fmax)) - ...
         EF(gf,df,S1,fmax,Ep(t),delta_t) - ...
         RF(td,S1));
     
    tmpf_S2 = @(S1,S2,S3) ...
        (QSSE(PHISS(ass,bss,S2,amin,amax),PHIC(ac,bc,S2,amin,amax),PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)) + ...
         QSIE(PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c),FA(max(0,PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c)*min(1,(1-PHISS(ass,bss,S2,amin,amax))/(1-PHIC(ac,bc,S2,amin,amax)))),FSS(af,bf,S3,bmax,S1,fmax))+QSSE(PHISS(ass,bss,S2,amin,amax),PHIC(ac,bc,S2,amin,amax),PC(PG(bg,ag,P(t))-QSE(ac,bc,S2,amin,amax,PG(bg,ag,P(t))),c))) + ...
         QB(bb,ab,S3,bmax) - ...
         EA1(S2,PHIC(ac,bc,S2,amin,amax)*Ep(t),delta_t) - ...
         EA2(ga,da,S2,amax,Ep(t),delta_t) - ...
         RA(PHIC(ac,bc,S2,amin,amax),FSS(af,bf,S3,bmax,S1,fmax),delta_t) - ...
         QA(aa,ba,S2,amin,amax,1));
     
    tmpf_S3 = @(S1,S2,S3) ...
        (RF(td,S1) + ...
         RA(PHIC(ac,bc,S2,amin,amax),FSS(af,bf,S3,bmax,S1,fmax),delta_t) - ...
         EB(gb,db,S3,bmax,Ep(t),delta_t) - ...
         QB(bb,ab,S3,bmax));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3);        % this returns a new anonymous function that we solve in the next step

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
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    tmp_phiss(t)= PHISS(ass,bss,tmp_sFlux(2),amin,amax);
    tmp_phic(t) = PHIC(ac,bc,tmp_sFlux(2),amin,amax);
    tmp_fss(t)  = FSS(af,bf,tmp_sFlux(3),bmax,tmp_sFlux(1),fmax);
    
    flux_pg(t)  = PG(bg,ag,P(t));
    flux_ei(t)  = EI(P(t),flux_pg(t));
    flux_qse(t) = QSE(ac,bc,tmp_sFlux(2),amin,amax,flux_pg(t));
    flux_pc(t)  = PC(flux_pg(t)-flux_qse(t),c);
    flux_qie(t) = QIE(flux_pg(t)-flux_qse(t),flux_pc(t));
    flux_qsse(t)= QSSE(tmp_phiss(t),tmp_phic(t),flux_pc(t));
    flux_fa(t)  = FA(max(0,flux_pc(t)*min(1,(1-tmp_phiss(t))/(1-tmp_phic(t)))),tmp_fss(t));
    flux_qsie(t)= QSIE(flux_pc(t),flux_fa(t)+flux_qsse(t));
    flux_ef(t)  = EF(gf,df,tmp_sFlux(1),fmax,Ep(t),delta_t);
    flux_rf(t)  = RF(td,tmp_sFlux(1));
    flux_ea1(t) = EA1(tmp_sFlux(2),tmp_phic(t)*Ep(t),delta_t) ;
    flux_ea2(t) = EA2(ga,da,tmp_sFlux(2),amax,Ep(t),delta_t);
    flux_qa(t)  = QA(aa,ba,tmp_sFlux(2),amin,amax,1);
    flux_ra(t)  = RA(tmp_phic(t),tmp_fss(t),delta_t);
    flux_qb(t)  = QB(bb,ab,tmp_sFlux(3),bmax);
    flux_eb(t)  = EB(gb,db,tmp_sFlux(3),bmax,Ep(t),delta_t);
    
    % Update the stores
    store_S1(t) = S1old + (flux_fa(t) - flux_ef(t) - flux_rf(t)) * delta_t;
    store_S2(t) = S2old + (flux_qsse(t) + flux_qsie(t) + flux_qb(t) - flux_ea1(t) - ...
                           flux_ea2(t) - flux_ra(t) - flux_qa(t)) * delta_t;
    store_S3(t) = S3old + (flux_rf(t) + flux_ra(t) - flux_eb(t) - flux_qb(t)) * delta_t; 
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ei + flux_ef + flux_ea1 + flux_ea2 + flux_eb) * delta_t;
    fluxOutput.Q      = (flux_qse + flux_qie + flux_qa) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ei   = flux_ei * delta_t;
    fluxInternal.ef   = flux_ef * delta_t;
    fluxInternal.ea   = (flux_ea1+flux_ea2) * delta_t;
    fluxInternal.eb   = flux_eb * delta_t;
    fluxInternal.pg   = flux_pg * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qie  = flux_qie * delta_t;
    fluxInternal.pc   = flux_pc * delta_t;
    fluxInternal.qsse = flux_qsse * delta_t;
    fluxInternal.qsie = flux_qsie * delta_t;
    fluxInternal.fa   = flux_fa * delta_t;
    fluxInternal.rf   = flux_rf * delta_t;
    fluxInternal.ra   = flux_ra * delta_t;
    fluxInternal.qa   = flux_qa * delta_t;
    fluxInternal.qb   = flux_qb * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


 




