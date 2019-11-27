function [ fluxOutput, fluxInternal, storeInternal, waterBalance  ] = ...
            m_14_topmodel_7p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: TOPMODEL
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Beven, K., Lamb, R., Quinn, P., Romanowicz, R., & Freer, J. (1995). 
% TOPMODEL. In V. P. Singh (Ed.), Computer Models of Watershed Hydrology 
% (pp. 627–668). Baton Rouge: Water Resources Publications, USA.
%
% Beven, K. J., & Kirkby, M. J. (1979). A physically based, variable 
% contributing area model of basin hydrology / Un modèle à base physique 
% de zone d’appel variable de l'hydrologie du bassin versant. Hydrological 
% Sciences Bulletin, 24(1), 43–69. http://doi.org/10.1080/02626667909491834
%
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. a., Vrugt, J. a., 
% Gupta, H. V., … Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between 
% hydrological models. Water Resources Research, 44(12). 
% http://doi.org/10.1029/2007WR006735
%
% Sivapalan, M., Beven, K., & Wood, E. F. (1987). On hydrologic similarity:
% 2. A scaled model of storm runoff production. Water Resources Research, 
% 23(12), 2266–2278. http://doi.org/10.1029/WR023i012p02266

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
suzmax  = theta(1);     % Maximum soil moisture storage in unsatured zone [mm]
st      = theta(2);     % Threshold for flow generation and evap change as fraction of suzmax [-]
kd      = theta(3);     % Leakage to saturated zone flow coefficient [mm/d]
q0      = theta(4);     % Zero deficit base flow speed [mm/d]
f       = theta(5);     % Baseflow scaling coefficient [mm-1]
chi     = theta(6);     % Gamma distribution parameter [-]
phi     = theta(7);     % Gamma distribution parameter [-]
mu      = 3;            % Gamma distribution parameter, fixed (Clark et al, 2008)
lambda  = chi*phi+mu;   % Ac computation parameter, mean of the gamma distribution

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial upper zone storage
S20     = storeInitial(2);       % Initial saturated zone storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];                   % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_qof  = zeros(1,t_end);
flux_peff = zeros(1,t_end);
flux_ea   = zeros(1,t_end);
flux_qex  = zeros(1,t_end);
flux_qv   = zeros(1,t_end);
flux_qb   = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Upper zone
% S2. Saturated zone

% QOF(chi,phi,3,lambda,f,S2,P(t)): direct runoff from saturated part
QOF = saturation_7;

% QEX(P(t)-QOF(chi,phi,3,lambda,f,S2,P(t)),S1,suzmax): unsaturated zone
% overflow
QEX = saturation_1;

% EA(st,S1,suzmax,Ep(t),delta_t); evaporation from unsaturated zone
EA = evap_3;

% QV(S1,kd,st*suzmax,suzmax-st*suzmax): threshold drainage to saturated zone
QV = interflow_10;

% QB(q0,f,S2): exponentially decaying flow from deficit store
QB = baseflow_4;

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
                        QOF(chi,phi,3,lambda,f,S2,P(t)) - ...
                        QEX(P(t)-QOF(chi,phi,3,lambda,f,S2,P(t)),S1,suzmax) - ...
                        EA(st,S1,suzmax,Ep(t),delta_t) - ...
                        QV(S1,kd,st*suzmax,suzmax-st*suzmax));              % store 1 function with current flux values
    tmpf_S2 = @(S1,S2) (QB(q0,f,S2) - ...
                        QV(S1,kd,st*suzmax,suzmax-st*suzmax));              % store 2 function (DEFICIT!)
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(eq_sys(1),eq_sys(2)),...    % system of storage equations
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
    flux_qof(t)  = QOF(chi,phi,3,lambda,f,tmp_sFlux(2),P(t));
    flux_peff(t) = P(t) - flux_qof(t);
    flux_ea(t)   = EA(st,tmp_sFlux(1),suzmax,Ep(t),delta_t);
    flux_qex(t)  = QEX(flux_peff(t),tmp_sFlux(1),suzmax);
    flux_qv(t)   = QV(tmp_sFlux(1),kd,st*suzmax,suzmax-st*suzmax);
    flux_qb(t)   = QB(q0,f,tmp_sFlux(2));
    
    % Update the stores
    store_S1(t) = S1old + (flux_peff(t) - flux_ea(t) - flux_qex(t) - flux_qv(t)) * delta_t;
    store_S2(t) = S2old + (flux_qb(t) - flux_qv(t)) * delta_t;    
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = (flux_qof + flux_qex + flux_qb) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.qof  = flux_qof * delta_t;
    fluxInternal.peff = flux_peff * delta_t;
    fluxInternal.ea   = flux_ea * delta_t;
    fluxInternal.qex  = flux_qex * delta_t;
    fluxInternal.qv   = flux_qv * delta_t;
    fluxInternal.qb   = flux_qb * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4
    
    waterBalance =  sum(P).*delta_t - ...           % Incoming precipitation
                    sum(fluxOutput.Q) - ...         % Outgoing flow
                    sum(fluxOutput.Ea) - ...        % Outgoing evaporation
                    (store_S1(end)-S10) + ...       % Store change
                    (store_S2(end)-S20);            % Store change (deficit store)

    disp(['Total P  = ',num2str(sum(P)),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm (deficit).'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(routing)) + delta S = ',num2str(waterBalance)]); 
end


 




