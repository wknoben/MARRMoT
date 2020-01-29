function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
        m_10_susannah2_6p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Susannah Brook model v2
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Son, K., & Sivapalan, M. (2007). Improving model structure and reducing 
% parameter uncertainty in conceptual water balance models through the use 
% of auxiliary data. Water Resources Research, 43(1). 
% https://doi.org/10.1029/2006WR005032

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
sb    = theta(1);     % Maximum soil moisture storage [mm]
phi   = theta(2);     % Porosity [-]
fc    = theta(3);     % Field capacity as fraction of sb [-]
r     = theta(4);     % Fraction of recharge coefficient [-]
c     = theta(5);     % Subsurface flow constant [1/d]
d     = theta(6);     % Subsurface flow constant [-]

%%INITIALISE MODEL STORES AND ROUTING VECTORS
S10   = storeInitial(1);       % Initial unsaturated zone storage
S20   = storeInitial(2);       % Initial saturated zone storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];             % lower bounds of stores
store_upp = [];                % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_eus  = zeros(1,t_end);
flux_rg   = zeros(1,t_end);
flux_se   = zeros(1,t_end);
flux_esat = zeros(1,t_end);
flux_qse  = zeros(1,t_end);
flux_qss  = zeros(1,t_end);
flux_qr   = zeros(1,t_end);
flux_qt   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Unsaturated store (variable capacity)
% S2. Saturated store

% EUS(S1,sb,Ep,delta_t): evaporation from unsaturated storage
EUS = evap_7;

% RG(P,S1,(sb-S2)*fc/phi): overflow to saturated zone
RG = saturation_1;

% SE(S1,(sb-S2)*fc/phi,delta_t): excess flow due to UZ storage change
SE = excess_1;

% ESAT(S2,sb,Ep,delta_t): evaporation from saturated storage
ESAT = evap_7; 

% QSE(RG+SE,S2,sb): saturation excess runoff
QSE = saturation_1;

% QSS((1-r)*c,d,S2,delta_t): non-linear subsurface flow
QSS = interflow_3;

% QR(r*c,d,S2,delta_t): recharge to deep groundwater (sink term)
QR = interflow_3;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,1;
                                               1,1]);                       % Specify the Jacobian pattern                                               
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
    tmpf_S1 = @(S1,S2) (P(t) - ...
                        EUS(S1,sb,Ep(t),delta_t) - ...
                        RG(P(t),S1,(sb-S2)*fc/phi) - ...
                        SE(S1,(sb-S2)*fc/phi,delta_t));                     % store 1 function with current flux values
    tmpf_S2 = @(S1,S2) (RG(P(t),S1,(sb-S2)*fc/phi) + ...
                        SE(S1,(sb-S2)*fc/phi,delta_t) - ...
                        ESAT(S2,sb,Ep(t),delta_t) - ...
                        QSE(RG(P(t),S1,(sb-S2)*fc/phi)+SE(S1,(sb-S2)*fc/phi,delta_t),S2,sb) - ...
                        QSS((1-r)*c,d,S2,delta_t) - ...
                        QR(r*c,d,S2,delta_t));                              % store 2 function
    
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
    flux_eus(t)  = EUS(tmp_sFlux(1),sb,Ep(t),delta_t);
    flux_rg(t)   = RG(P(t),tmp_sFlux(1),(sb-tmp_sFlux(2))*fc/phi);
    flux_se(t)   = SE(tmp_sFlux(1),(sb-tmp_sFlux(2))*fc/phi,delta_t);
    flux_esat(t) = ESAT(tmp_sFlux(2),sb,Ep(t),delta_t);
    flux_qse(t)  = QSE(flux_rg(t)+flux_se(t),tmp_sFlux(2),sb);
    flux_qss(t)  = QSS((1-r)*c,d,tmp_sFlux(2),delta_t);
    flux_qr(t)   = QR(r*c,d,tmp_sFlux(2),delta_t);
    flux_qt(t)   = flux_qse(t) + flux_qss(t);
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_eus(t) - flux_rg(t) - flux_se(t)) * delta_t;
    store_S2(t) = S2old + (flux_rg(t) + flux_se(t) - flux_esat(t) - ...
                    flux_qse(t) - flux_qss(t) - flux_qr(t)) * delta_t;    
   
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_eus + flux_esat) * delta_t;
    fluxOutput.Qr     = flux_qr * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.eus  = flux_eus * delta_t;
    fluxInternal.rg   = flux_rg * delta_t;
    fluxInternal.se   = flux_se * delta_t;
    fluxInternal.esat = flux_esat * delta_t;
    fluxInternal.qse  = flux_qse * delta_t;
    fluxInternal.qss  = flux_qss * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    waterBalance =  sum(P).*delta_t - ...                           % Incoming precipitation 
                    sum(fluxOutput.Q) - ...                         % Flux Q leaving the model
                    sum(fluxOutput.Ea) - ...                        % Flux Ea leaving the model
                    sum(fluxOutput.Qr) - ...                        % Groundwater sink term
                    ((store_S1(end)-S10) + (store_S2(end)-S20));    % Storage change 
        
    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Total Q_sink  = ',num2str(sum(fluxOutput.Qr)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(Qg) + sum(E_a) + delta S) = ',num2str(waterBalance)]);                          % delta S

end



 




