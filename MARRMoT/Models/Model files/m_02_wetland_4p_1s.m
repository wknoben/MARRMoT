function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_02_wetland_4p_1s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Wetland model 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Savenije, H. H. G. (2010). “Topography driven conceptual modelling 
% (FLEX-Topo).” Hydrology and Earth System Sciences, 14(12), 2681–2692. 
% https://doi.org/10.5194/hess-14-2681-2010


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
t_end = length(P);

% Parameters
% [name in documentation] = theta(order in which specified in parameter file)
dw    = theta(1);     % Daily interception [mm]
betaw = theta(2);     % Soil moisture storage distribution parameter [-]
swmax = theta(3);     % Maximum soil moisture storage [mm], 
kw    = theta(4);     % Runoff coefficient [d-1]

%%INITIALISE MODEL STORES AND ROUTING VECTORS
S10   = storeInitial(1);            % Initial soil moisture storage

%%DEFINE STORE BOUNDARIES
store_min = [0];                    % lower bounds of stores
store_upp = [];                     % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);

flux_pe    = zeros(1,t_end);
flux_ei    = zeros(1,t_end);
flux_ew    = zeros(1,t_end);
flux_qwsof = zeros(1,t_end);
flux_qwgw  = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% With Matlab's fsolver, smoothing is only needed when the function is
% undefined (i.e. thresholds). Angle discontinuities (such as from
% min(0,x)) can be dealt with by the solver.

% Store numbering:
% S1. Soil moisture

% PE(P(t),dw): effective rainfall after interception
PE = interception_2;

% EW(S1,Ep,delta_t): evaporation from soil moisture
EW = evap_1;

% QWSOF(S1,swmax,betaw,PE(P(t),dw)): saturation excess flow
QWSOF = saturation_2;

% QWGW(kw,S1): baseflow
QWGW = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fzero_options = optimset('Display','off');                                  % settings of the root finding method
lsqnonlin_options = optimoptions('lsqnonlin','Display','none',...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1)    (PE(P(t),dw) - ...
                        EW(S1,Ep(t),delta_t) - ...
                        QWSOF(S1,swmax,betaw,PE(P(t),dw)) - ...
                        QWGW(kw,S1));                                       % store 1 function with current flux values
        
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old],...
                      delta_t,...
                      tmpf_S1);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fzero(solve_fun,S1old,fzero_options);

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(eq_sys(1)),...  % system of ODEs
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old], ...                        % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
      
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_pe(t)    = PE(P(t),dw);
    flux_ei(t)    = P(t) - flux_pe(t);
    flux_ew(t)    = EW(tmp_sFlux(1),Ep(t),delta_t);
    flux_qwsof(t) = QWSOF(tmp_sFlux(1),swmax,betaw,flux_pe(t));
    flux_qwgw(t)  = QWGW(kw,tmp_sFlux(1));
    
    % Update the stores
    store_S1(t) = S1old + (flux_pe(t) - flux_ew(t) - flux_qwsof(t) - flux_qwgw(t)) * delta_t;
       
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    fluxOutput.Ea     = (flux_ew + flux_ei) * delta_t;
    fluxOutput.Q      = (flux_qwsof + flux_qwgw) * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pe     = flux_pe * delta_t;
    fluxInternal.ei     = flux_ei * delta_t;
    fluxInternal.ew     = flux_ew * delta_t;
    fluxInternal.Qwsof  = flux_qwsof * delta_t;
    fluxInternal.Qwgw   = flux_qwgw * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;

% Check water balance
if nargout == 4
    waterBalance = checkWaterBalance(P.*delta_t,...     % Incoming precipitation
                                     fluxOutput,...     % Fluxes Q and Ea leaving the model
                                     storeInternal,...  % Time series of storages ...
                                     storeInitial,...   % And initial store values to calculate delta S
                                     0);                % Whether the model uses a routing scheme that
                                                        % still contains water. Use '0' for no routing
end


end
 




