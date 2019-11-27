function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_46_classic_12p_8s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: [CLASSIC] 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Crooks, S. M., & Naden, P. S. (2007). CLASSIC: a semi-distributed 
% rainfall-runoff modelling system. Hydrology and Earth System Sciences, 
% 11(1), 516–531. http://doi.org/10.5194/hess-11-516-2007

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
fap   = theta(1);     % Fraction of catchment area that has permeable soils [-]
fdp   = theta(2);     % Fraction of depth of permeable soil that is store Px [-]
dp    = theta(3);     % Depth of permeable soil [mm]
cq    = theta(4);     % Runoff coefficient for permeable soil [d-1]
d1    = theta(5);     % Fraction of Ps that infiltrates into semi-permeable soil [-]
tf    = theta(6);     % Fraction of (1-fap) that is fas [-]
fds   = theta(7);     % Fraction of depth of semi-permeable soil that is store Sx [-]
ds    = theta(8);     % Depth of semi-permeable soil [mm]
d2    = theta(9);     % Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]
cxq   = theta(10);    % Quick runoff coefficient for semi-permeable soil [d-1]
cxs   = theta(11);    % Slow runoff coefficient for semi-permeable soil [d-1]
cu    = theta(12);    % Runoff coefficient for impermeable soil [d-1]

% Auxiliary parameters
fas = (1-fap)*tf;     % Fraction of catchment area that has semi-permeable soils [-]
fai = 1-fap-fas;      % Fraction of catchment area that has impermeable soils [-]
pxm = fdp*dp;         % Depth of store Px [mm]
pym = (1-fdp)*dp;     % Depth of store Py [mm]
sxm = fds*ds;         % Depth of store Sx [mm]
sym = (1-fds)*ds;     % Depth of store Sy [mm]

%%INITIALISE MODEL STORES
S10 = storeInitial(1);       % Initial permeable Px storage
S20 = storeInitial(2);       % Initial permeable Py storage
S30 = storeInitial(3);       % Initial permeable routing
S40 = storeInitial(4);       % Initial semi-permeable Sx storage
S50 = storeInitial(5);       % Initial semi-permeable Sy storage
S60 = storeInitial(6);       % Initial semi-permeable quick routing
S70 = storeInitial(7);       % Initial semi-permeable slow routing
S80 = storeInitial(8);       % Initial impermeable routing

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0,0,0,0,0];       % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS (all upper case)
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);
store_S5 = zeros(1,t_end);
store_S6 = zeros(1,t_end);
store_S7 = zeros(1,t_end);
store_S8 = zeros(1,t_end);

flux_pp  = zeros(1,t_end);
flux_ps  = zeros(1,t_end);
flux_pi  = zeros(1,t_end);
flux_epx = zeros(1,t_end);
flux_ppx = zeros(1,t_end);
flux_epy = zeros(1,t_end);
flux_ppe = zeros(1,t_end);
flux_q   = zeros(1,t_end);
flux_psd = zeros(1,t_end);
flux_psi = zeros(1,t_end);
flux_esx = zeros(1,t_end);
flux_psx = zeros(1,t_end);
flux_esy = zeros(1,t_end);
flux_pse = zeros(1,t_end);
flux_psq = zeros(1,t_end);
flux_pss = zeros(1,t_end);
flux_xq  = zeros(1,t_end);
flux_xs  = zeros(1,t_end);
flux_ei  = zeros(1,t_end);
flux_pie = zeros(1,t_end);
flux_u   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% With Matlab's fsolver, smoothing is only needed when the function is
% undefined (i.e. thresholds). Angle discontinuities (such as from
% min(0,x)) can be dealt with by the solver.

% Store numbering:
% S1. Permeable soil Px
% S2. Permeable soil deficit Py
% S3. Permeable soil routing
% S4. Semi-permeable soil Px
% S5. Semi-permeable soil deficit Py
% S6. Semi-permeable soil quick routing
% S7. Semi-permeable soil slow routing
% S8. Impermeable soil routing

% PP(fap,P(t)): fraction P on permeable soil
PP = split_1;

% PS(fas,P(t)): fraction P on semi-permeable soil
PS = split_1;

% PI(fai,P(t)): fraction P on impermeable soil
PI = split_1;

% EPX(S1,fap*Ep(t),delta_t): evaporation from Px (S1)
EPX = evap_1;

% PPX(PP(fap,P(t)),S1,pxm): overflow from Px
PPX = saturation_1;

% EPY(1.9,0.6523,pxm,S2+pxm,fap*Ep(t)-EPX(S1,fap*Ep(t),delta_t)): evaporation from Py
EPY = evap_18;

% PPE(PPX(PP(fap,P(t)),S1,pxm),S2,0.01): overflow from Py
PPE = saturation_9;

% Q(cq,S3): baseflow from permeable soils
Q = baseflow_1;

% PSI(d1,PS(fas,P(t))): fraction of PS that is infiltration
PSI = split_1;

% PSD(1-d1,PS(fas,P(t))): fraction of PS that is direct runoff
PSD = split_1;

% ESX(S4,fas*Ep(t),delta_t): evaporation from Sx (S4)
ESX = evap_1;

% PSX(PSI(d1,PS(fas,P(t))),S4,sxm): overflow from Sx
PSX = saturation_1;

% ESY(1.9,0.6523,sxm,S4+sxm,fas*Ep(t)-ESX(S4,fas*Ep(t))): evaporation from Sy
ESY = evap_18;

% PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,0.01): overflow from Sy
PSE = saturation_9;

% PSQ(d2,PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,sym)+PSD(1-d1,PS(fas,P(t)))):
% effective flow from semi-permeable area into quick flow
PSQ = split_1;

% PSS(1-d2,PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,sym)+PSD(1-d1,PS(fas,P(t)))):
% effective flow from semi-permeable area into slow flow
PSS = split_1;

% XQ(cxq,S6): quick flow from semi-permeable soils
XQ = baseflow_1;

% XS(cxs,S7): slow flow from semi-permeable soils
XS = baseflow_1;

% PIE(fai*PI(fai,P(t)),0.5):  effective precipitation on impermeable soil
PIE = effective_1;

% U(cu,S8): flow from impermeable soils
U = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0,0,0,0,0;
                                               1,1,0,0,0,0,0,0;
                                               1,1,1,0,0,0,0,0;
                                               0,0,0,1,0,0,0,0;
                                               0,0,0,1,1,0,0,0;
                                               0,0,0,1,1,1,0,0;
                                               0,0,0,1,1,0,1,0;
                                               0,0,0,0,0,0,0,1]);           % Specify the Jacobian pattern                                               
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0,0,0,0,0;
                                                  1,1,0,0,0,0,0,0;
                                                  1,1,1,0,0,0,0,0;
                                                  0,0,0,1,0,0,0,0;
                                                  0,0,0,1,1,0,0,0;
                                                  0,0,0,1,1,1,0,0;
                                                  0,0,0,1,1,0,1,0;
                                                  0,0,0,0,0,0,0,1],...
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
    if t == 1; S8old = S80; else; S8old = store_S8(t-1); end                % store 8 at t-1

    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PP(fap,P(t)) - ...
                EPX(S1,fap*Ep(t),delta_t) - ...
                PPX(PP(fap,P(t)),S1,pxm));
            
    tmpf_S2 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...                                % Deficit store!
               (-PPX(PP(fap,P(t)),S1,pxm) + ...
                EPY(1.9,0.6523,pxm,S2+pxm,fap*Ep(t)-EPX(S1,fap*Ep(t),delta_t)) + ...
                PPE(PPX(PP(fap,P(t)),S1,pxm),S2,0.01));
    
    tmpf_S3 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PPE(PPX(PP(fap,P(t)),S1,pxm),S2,0.01) - ...
                Q(cq,S3));
    
    tmpf_S4 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PSI(d1,PS(fas,P(t))) - ...
                ESX(S4,fas*Ep(t),delta_t) - ...
                PSX(PSI(d1,PS(fas,P(t))),S4,sxm));
            
    tmpf_S5 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...                                % Deficit store!
               (-PSX(PSI(d1,PS(fas,P(t))),S4,sxm) + ...
                ESY(1.9,0.6523,sxm,S4+sxm,fas*Ep(t)-ESX(S4,fas*Ep(t),delta_t)) + ...
                PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,0.01));
            
    tmpf_S6 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PSQ(d2,PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,0.01)+PSD(1-d1,PS(fas,P(t)))) - ...
                XQ(cxq,S6));

    tmpf_S7 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PSS(1-d2,PSE(PSX(PSI(d1,PS(fas,P(t))),S4,sxm),S5,0.01)+PSD(1-d1,PS(fas,P(t)))) - ...
                XS(cxs,S7));
            
    tmpf_S8 = @(S1,S2,S3,S4,S5,S6,S7,S8) ...
               (PIE(fai*PI(fai,P(t)),0.5) - ...
                U(cu,S8));
            
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old,S4old,...
                       S5old,S6old,S7old,S8old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4,...
                      tmpf_S5,tmpf_S6,tmpf_S7,tmpf_S8);                     % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                        eq_sys(5),eq_sys(6),eq_sys(7),eq_sys(8)),...        % system of storage equations
                        [S1old,S2old,S3old,S4old,...
                         S5old,S6old,S7old,S8old],...                       % storage values on previous time step
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
                                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4),...
                                        eq_sys(5),eq_sys(6),eq_sys(7),eq_sys(8)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old,...
                                         S5old,S6old,S7old,S8old], ...      % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end

% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_pp(t)  = PP(fap,P(t));
    flux_ps(t)  = PS(fas,P(t));
    flux_pi(t)  = PI(fai,P(t));
    flux_epx(t) = EPX(tmp_sFlux(1),fap*Ep(t),delta_t);
    flux_ppx(t) = PPX(flux_pp(t),tmp_sFlux(1),pxm);
    flux_epy(t) = EPY(1.9,0.6523,pxm,tmp_sFlux(2)+pxm,fap*Ep(t)-flux_epx(t));
    flux_ppe(t) = PPE(flux_ppx(t),tmp_sFlux(2),0.01);
    flux_q(t)   = Q(cq,tmp_sFlux(3));
    flux_psd(t) = PSD(1-d1,flux_ps(t));
    flux_psi(t) = PSI(d1,flux_ps(t));
    flux_esx(t) = ESX(tmp_sFlux(4),fas*Ep(t),delta_t);
    flux_psx(t) = PSX(flux_psi(t),tmp_sFlux(4),sxm);
    flux_esy(t) = ESY(1.9,0.6523,sxm,tmp_sFlux(4)+sxm,fas*Ep(t)-flux_esx(t));
    flux_pse(t) = PSE(flux_psx(t),tmp_sFlux(5),0.01);
    flux_psq(t) = PSQ(d2,flux_pse(t)+flux_psd(t));
    flux_pss(t) = PSS(1-d2,flux_pse(t)+flux_psd(t));
    flux_xq(t)  = XQ(cxq,tmp_sFlux(6));
    flux_xs(t)  = XS(cxs,tmp_sFlux(7));
    flux_pie(t) = PIE(fai*flux_pi(t),0.5);
    flux_ei(t)  = flux_pi(t) - flux_pie(t);
    flux_u(t)   = U(cu,tmp_sFlux(8));
    
    % Update the stores
    store_S1(t) = S1old + (flux_pp(t)  - flux_epx(t) - flux_ppx(t)) * delta_t;
    store_S2(t) = S2old - (flux_ppx(t) - flux_epy(t) - flux_ppe(t)) * delta_t;    
    store_S3(t) = S3old + (flux_ppe(t) - flux_q(t)) * delta_t;    
    store_S4(t) = S4old + (flux_psi(t) - flux_esx(t) - flux_psx(t)) * delta_t;    
    store_S5(t) = S5old - (flux_psx(t) - flux_esy(t) - flux_pse(t)) * delta_t;    
    store_S6(t) = S6old + (flux_psq(t) - flux_xq(t)) * delta_t;    
    store_S7(t) = S7old + (flux_pss(t) - flux_xs(t)) * delta_t;    
    store_S8(t) = S8old + (flux_pie(t) - flux_u(t)) * delta_t;    
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_epx + flux_epy + flux_esx + flux_esy + flux_ei) * delta_t;
    fluxOutput.Q      = (flux_q + flux_xq + flux_xs + flux_u) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pp  = flux_pp  * delta_t;
    fluxInternal.ps  = flux_ps  * delta_t;
    fluxInternal.pi  = flux_pi  * delta_t;
    fluxInternal.epx = flux_epx * delta_t;
    fluxInternal.ppx = flux_ppx * delta_t;
    fluxInternal.epy = flux_epy * delta_t;
    fluxInternal.ppe = flux_ppe * delta_t;
    fluxInternal.q   = flux_q   * delta_t;
    fluxInternal.psd = flux_psd * delta_t;
    fluxInternal.esx = flux_esx * delta_t;
    fluxInternal.psx = flux_psx * delta_t;
    fluxInternal.esy = flux_esy * delta_t;
    fluxInternal.pse = flux_pse * delta_t;
    fluxInternal.psq = flux_psq * delta_t;
    fluxInternal.pss = flux_pss * delta_t;
    fluxInternal.xq  = flux_xq  * delta_t;
    fluxInternal.xs  = flux_xs  * delta_t;
    fluxInternal.pie = flux_pie * delta_t;
    fluxInternal.ei  = flux_ei  * delta_t;
    fluxInternal.u   = flux_u   * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    storeInternal.S5  = store_S5;
    storeInternal.S6  = store_S6;
    storeInternal.S7  = store_S7;
    storeInternal.S8  = store_S8;

% Check water balance
if nargout == 4
    
    % Manual water balance on account of deficit stores
    waterBalance = ...
        sum(P).*delta_t - ...           % Incoming precipitation
        sum(fluxOutput.Ea) - ...        % Outgoing evaporation
        sum(fluxOutput.Q) - ...         % Outgoing flow
        (store_S1(end)-S10) + ...       % Storage change
        (store_S2(end)-S20) - ...           % Deficit store
        (store_S3(end)-S30) - ...
        (store_S4(end)-S40) + ...
        (store_S5(end)-S50) - ...           % Deficit store
        (store_S6(end)-S60) - ...
        (store_S7(end)-S70) - ...
        (store_S8(end)-S80);

    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm deficit.'])
    disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
    disp(['Delta S4 = ',num2str((store_S4(end)-S40)),'mm.'])
    disp(['Delta S5 = ',num2str((store_S5(end)-S50)),'mm deficit.'])
    disp(['Delta S6 = ',num2str((store_S6(end)-S60)),'mm.'])
    disp(['Delta S7 = ',num2str((store_S7(end)-S70)),'mm.'])
    disp(['Delta S8 = ',num2str((store_S8(end)-S80)),'mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a)) - delta S = ',...
        num2str(waterBalance)]);

end
 




