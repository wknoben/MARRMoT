function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_05_ihacres_7p_1s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: IHACRES
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Croke, B. F. W., & Jakeman, A. J. (2004). A catchment moisture deficit 
% module for the IHACRES rainfall-runoff model. Environmental Modelling and
% Software, 19(1), 1–5. http://doi.org/10.1016/j.envsoft.2003.09.001
%
% Littlewood, I. G., Down, K., Parker, J. R., & Post, D. A. (1997). IHACRES
% v1.0 User Guide.
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral 
% catchments. Water Resources Research, 33(1), 153–166. 
% http://doi.org/doi:10.1029/96WR02840


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
lp      = theta(1);     % Wilting point [mm]
d       = theta(2);     % Threshold for flow generation [mm]
p       = theta(3);     % Flow response non-linearity [-]
alpha   = theta(4);     % Fast/slow flow division [-]
tau_q   = theta(5);     % Fast flow routing delay [d]
tau_s   = theta(6);     % Slow flow routing delay [d]
tau_d   = theta(7);     % Pure time delay of total flow [d]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial soil moisture deficit

%%DEFINE STORE BOUNDARIES
store_min = [0];                     % lower bounds of stores
store_upp = [];                      % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);

flux_ea   = zeros(1,t_end);
flux_u    = zeros(1,t_end);
flux_uq   = zeros(1,t_end);
flux_us   = zeros(1,t_end);
flux_xq   = zeros(1,t_end);
flux_xs   = zeros(1,t_end);
flux_qt   = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_q] = uh_5_half(1,tau_q,delta_t);
[~,uh_s] = uh_5_half(1,tau_s,delta_t);
[~,uh_t] = uh_8_delay(1,tau_d,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_xq_old  = zeros(1,length(uh_q));  % temporary vector needed to deal with routing
tmp_xs_old  = zeros(1,length(uh_s));  % temporary vector needed to deal with routing
tmp_qt_old  = zeros(1,length(uh_t));  % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Current moisture deficit store

% EA(S1,lp,Ep(t)): evaporation from soil moisture
EA = evap_12;

% U(S1,d,p,P(t)): effective precipitation
U = saturation_5;

% UQ(alpha,U(S1,d,p,P(t))): part of effective precipitation to fast flow
UQ = split_1;

% US(1-alpha,U(S1,d,p,P(t))): part of effective precipitation to fast flow
US = split_1;

%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme  = solver.name;                                                      

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fzero_options = optimset('Display','off');                                  % [1 store] settings of the root finding method
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'MaxFunEvals',1000);

%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = @(S1) (-P(t) + ...
                     EA(S1,lp,Ep(t)) +...
                     U(S1,d,p,P(t)));                                       % store 1 function with current flux values
    
    % Call the numerical scheme function to create the ODE approximations
    solve_fun = feval(scheme,...
                      [S1old],...
                      delta_t,...
                      tmpf_S1);        % this returns a new anonymous function that we solve in the next step

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew, tmp_fval] = fzero(solve_fun,...
                                    S1old,...
                                    fzero_options);                         % 1 store solver

    % --- Check if the solver has found an acceptable solution and re-run
    % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
    % more robust. It runs solver.resnorm_iterations times, with different
    % starting points for the solver on each iteration ---
    tmp_resnorm = sum(tmp_fval.^2);
     
    if tmp_resnorm > solver.resnorm_tolerance
        [tmp_sNew,~,~] = rerunSolver('lsqnonlin', ...                       % [tmp_sNew,tmp_resnorm,flag]
                                        lsqnonlin_options, ...              % solver options
                                        @(eq_sys) solve_fun(...             % system of ODEs
                                                    eq_sys(1)), ...
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
    flux_ea(t)   = EA(tmp_sFlux(1),lp,Ep(t));
    flux_u(t)    = U(tmp_sFlux(1),d,p,P(t));
    flux_uq(t)   = UQ(alpha,flux_u(t));
    flux_us(t)   = US(1-alpha,flux_u(t));
    
    % Update the stores
    store_S1(t) = S1old - (P(t) - flux_ea(t) - flux_u(t)).*delta_t;

    % Apply routing of fast and slow flow
    tmp_xq_cur      = flux_uq(t).*uh_q;                                     % find how the runoff of this time step will be distributed in time
    tmp_xq_old      = tmp_xq_old + tmp_xq_cur;                              % update the 'still-to-flow-out' vector
    flux_xq(t)      = tmp_xq_old(1);                                        % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_xq_old      = circshift(tmp_xq_old,-1);                             % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_xq_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_qt(t) and is now shifted to the end of the 'still-to-flow-out' vector)

    tmp_xs_cur      = flux_us(t).*uh_s;                                     % find how the runoff of this time step will be distributed in time
    tmp_xs_old      = tmp_xs_old + tmp_xs_cur;                              % update the 'still-to-flow-out' vector
    flux_xs(t)      = tmp_xs_old(1);                                        % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_xs_old      = circshift(tmp_xs_old,-1);                             % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_xs_old(end) = 0;
    
    % Apply time delay
    tmp_qt_cur      = (flux_xq(t)+flux_xs(t)).*uh_t;
    tmp_qt_old      = tmp_qt_old + tmp_qt_cur;
    flux_qt(t)      = tmp_qt_old(1);
    tmp_qt_old      = circshift(tmp_qt_old,-1);
    tmp_qt_old(end) = 0;
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = flux_ea * delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.u    = flux_u * delta_t;
    fluxInternal.uq   = flux_uq * delta_t;
    fluxInternal.us   = flux_us * delta_t;
    fluxInternal.xq   = flux_xq * delta_t;
    fluxInternal.xs   = flux_xs * delta_t;
    
    % --- Stores ---
    storeInternal.S1  = store_S1;

% Check water balance
if nargout == 4
    
    % Water balance
    waterBalance = sum(P.*delta_t) - ...                                        % Incoming precipitation
                   sum(fluxOutput.Q) - ...                                      % Flux Q leaving the model
                   sum(fluxOutput.Ea) + ...                                     % Flux Ea leaving the model
                   (store_S1(end)-S10) - ...                                    % Deficit change
                   (sum(tmp_xq_old)+sum(tmp_xs_old)+sum(tmp_qt_old)).*delta_t;  % Water still in routing schemes

    % Display
    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])                             
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])                  
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])                 
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])                
    disp(['On route = ',num2str((sum(tmp_xq_old)+sum(tmp_xs_old)+sum(tmp_qt_old)).*delta_t),'mm.'])    
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(routing)) + delta S = ',num2str(waterBalance)]);   

end                                            


 




