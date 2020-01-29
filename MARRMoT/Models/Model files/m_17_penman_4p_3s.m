function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_17_penman_4p_3s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Penman
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Penman, H. L. (1950). the Dependence of Transpiration on Weather and Soil
% Conditions. Journal of Soil Science, 1(1), 74–89. 
% http://doi.org/10.1111/j.1365-2389.1950.tb00720.x
%
% Wagener, T., Lees, M. J., & Wheater, H. S. (2002). A toolkit for the 
% development and application of parsimonious hydrological models. In Singh,
% Frevert, & Meyer (Eds.), Mathematical Models of Small Watershed Hydrology
% - Volume 2 (pp. 91–139). Water Resources Publications LLC, USA.


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
smax  = theta(1);     % Maximum soil moisture storage [mm]
phi   = theta(2);     % Fraction of direct runoff [-]
gam   = theta(3);     % Evaporation reduction in lower zone [-]
k1    = theta(4);     % Runoff coefficient [d-1]

%%INITIALISE MODEL STORES
S10   = storeInitial(1);       % Initial root zone storage
S20   = storeInitial(2);       % Initial soil moisture deficit
S30   = storeInitial(3);       % Initial routing storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0];                 % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);

flux_ea   = zeros(1,t_end);
flux_qex  = zeros(1,t_end);
flux_u1   = zeros(1,t_end);
flux_q12  = zeros(1,t_end);
flux_et   = zeros(1,t_end);
flux_u2   = zeros(1,t_end);
flux_q    = zeros(1,t_end);

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Root zone
% S2. Soil moisture deficit
% S3. Routing reservoir

% EA(S1,Ep(t),delta_t): evap from root zone
EA = evap_1;

% QEX(P(t),S1,smax): overflow from soil moisture
QEX = saturation_1;

% U1(phi,QEX(P(t),S1,smax)): bypass towards routing reservoir
U1 = split_1;

% Q12(1-phi,QEX(P(t),S1,smax)): replenishing of deficit store
Q12 = split_1;

% ET(gam,S2,S1,0.01,Ep(t),delta_t): evaporation from deficit store
ET = evap_16;

% U2(Q12(1-phi,QEX(P(t),S1,smax)),S2,0.01): overflow from deficit store
U2 = saturation_9;

% Q(k1,S3): flow from routing store
Q = baseflow_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0;
                                               1,1,0;
                                               1,1,1]);                     % Specify the Jacobian pattern                                               

lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0;
                                                  1,1,0;
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
                       (P(t) - ...
                        EA(S1,Ep(t),delta_t) - ...
                        QEX(P(t),S1,smax));                                 

    tmpf_S2 = @(S1,S2,S3) ...
                        (ET(gam,S2,S1,0.01,Ep(t),delta_t) + ...
                         U2(Q12(1-phi,QEX(P(t),S1,smax)),S2,0.01) - ...
                         Q12(1-phi,QEX(P(t),S1,smax)));
    
    tmpf_S3 = @(S1,S2,S3) ...
                        (U1(phi,QEX(P(t),S1,smax)) + ...
                         U2(Q12(1-phi,QEX(P(t),S1,smax)),S2,0.01) - ...
                         Q(k1,S3));
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old,S2old,S3old],...
                      delta_t,...
                      tmpf_S1,tmpf_S2,tmpf_S3);                             % this returns a new anonymous function that we solve in the next step

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
    flux_ea(t)   = EA(tmp_sFlux(1),Ep(t),delta_t);
    flux_qex(t)  = QEX(P(t),tmp_sFlux(1),smax);
    flux_u1(t)   = U1(phi,flux_qex(t));
    flux_q12(t)  = Q12(1-phi,flux_qex(t));
    flux_et(t)   = ET(gam,tmp_sFlux(2),tmp_sFlux(1),0.01,Ep(t),delta_t);
    flux_u2(t)   = U2(flux_q12(t),tmp_sFlux(2),0.01);
    flux_q(t)    = Q(k1,tmp_sFlux(3));
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ea(t) - flux_qex(t)) * delta_t;
    store_S2(t) = S2old + (flux_et(t) + flux_u2(t) - flux_q12(t)) * delta_t;    
    store_S3(t) = S3old + (flux_u1(t) + flux_u2(t) - flux_q(t)) * delta_t;

end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_ea + flux_et) * delta_t;
    fluxOutput.Q      = flux_q * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ea   = flux_ea * delta_t;
    fluxInternal.qex  = flux_qex * delta_t;
    fluxInternal.u1   = flux_u1 * delta_t;
    fluxInternal.q12  = flux_q12 * delta_t;
    fluxInternal.et   = flux_et * delta_t;
    fluxInternal.u2   = flux_u2 * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;

% Check water balance
if nargout == 4
    
    waterBalance = sum(P).*delta_t-...          % Incoming precipitation                                                 
                   sum(fluxOutput.Ea) - ...     % Outgoing evaporation
                   sum(fluxOutput.Q) - ...      % Outgoing flow
                   (store_S1(end)-S10) + ...    % Storage change 
                   (store_S2(end)-S20) - ...    % Storage change (deficit store)
                   (store_S3(end)-S30);         % Storage change
    
    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 (deficit) = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(routing)) + delta S = ',...
            num2str(waterBalance)]);

end


    

    

 




