function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
            m_07_gr4j_4p_2s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: GR4J
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Perrin, C., Michel, C., & Andréassian, V. (2003). Improvement of a 
% parsimonious model for streamflow simulation. Journal of Hydrology, 
% 279(1-4), 275–289. http://doi.org/10.1016/S0022-1694(03)00225-7
%
% Santos, L., Thirel, G., & Perrin, C. (2017). State-space representation 
% of a bucket-type rainfall-runoff model: a case study with State-Space GR4
% (version 1.0). Geoscientific Model Development Discussions, 1–22. 
% http://doi.org/10.5194/gmd-2017-264


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
x1      = theta(1);     % Maximum soil moisture storage [mm]
x2      = theta(2);     % Water exchange coefficient [mm/d]
x3      = theta(3);     % Maximum routing store storage [mm]
x4      = theta(4);     % Flow delay [d]

%%INITIALISE MODEL STORES
S10     = storeInitial(1);       % Initial soil moisture storage
S20     = storeInitial(2);       % Initial routing storage

%%DEFINE STORE BOUNDARIES
store_min = [0,0];               % lower bounds of stores
store_upp = [x1,x3];             % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);

flux_pn   = zeros(1,t_end);
flux_ef   = zeros(1,t_end);
flux_en   = zeros(1,t_end);
flux_ps   = zeros(1,t_end);
flux_es   = zeros(1,t_end);
flux_perc = zeros(1,t_end);
flux_q9   = zeros(1,t_end);
flux_q1   = zeros(1,t_end);
flux_qr   = zeros(1,t_end);
flux_fr   = zeros(1,t_end);
flux_fq   = zeros(1,t_end);
flux_qt   = zeros(1,t_end);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_half] = uh_1_half(1,x4,delta_t);
[~,uh_full] = uh_2_full(1,2*x4,delta_t);

%%INITIALISE ROUTING VECTORS
tmp_Q9_old  = zeros(1,length(uh_half));  	% temporary vector needed to deal with routing
tmp_Q1_old  = zeros(1,length(uh_full));     % temporary vector needed to deal with routing

%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Soil moisture
% S2. Routing store

% PS(S1,x1,max(P(t)-Ep(t),0)): part of net precipitation to soil moisture.
PS = saturation_4;

% ES(S1,x1,max(Ep(t)-P(t),0)): evaporation based on net Ep
ES = evap_11;

% PERC(S1,x1): percolation to groundwater
PERC = percolation_3;

% F(3.5,S2,x3,x2): water exchange
F = recharge_2;

% QR(S2,x3): routed flow
QR = baseflow_3;

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
    
% Special note: this model has two lag functions that delay the flows into 
% the fast and slow reservoirs. Due to the current setup of the routing
% functions, store 1 (soil moisture) needs to be computed first. The new S1
% value is then used for both routing schemes to determine inflow into
% store 2 (fast reservoir) and store 3(slow reservoir) respectively. To
% keep the model file consistent with other model files we also need to
% initialise the 'tmp_sNew' vector:
tmp_sNew = NaN.*zeros(1,2);    

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1

    % Net climatology input
    flux_pn(t) = max(P(t)-Ep(t),0);
    flux_en(t) = max(Ep(t)-P(t),0);
    flux_ef(t) = P(t) - flux_pn(t);

    % --- Solve store 1 (soil moisture) -----------------------------------
        % Create temporary store ODE's that need to be solved
        tmpf_S1 = @(S1,S2) (PS(S1,x1,flux_pn(t)) - ...
                            ES(S1,x1,flux_en(t)) - ...
                            PERC(S1,x1));                                   % store 1 function with current flux values
                        
        % Call the numerical scheme function to create the ODE approximations
        solve_fun1 = feval(scheme,...
                        [S1old],...
                        delta_t,...
                        tmpf_S1);                
                        
        % Use the specified numerical scheme to solve storages
        [tmp_sNew(1), tmp_fval] = fzero(solve_fun1,...
                                    S1old,...
                                    fzero_options);                         % 1 store solver                 
                        
        % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(1),~,~] = rerunSolver('lsqnonlin', ...                % [tmp_sNew,tmp_resnorm,flag]
                                            lsqnonlin_options, ...          % solver options
                                            @(eq_sys) solve_fun1(...        % system of ODEs
                                                        eq_sys(1)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(1), ...                % recent estimates
                                            [S1old], ...                    % storages ate previous time step
                                            store_min(1), ...               % lower bounds
                                            store_upp(1));                  % upper bounds              
        end               
    
    % --- Handle the flow routing based on store 1 estimates --------------    
        % Find the storages needed to update fluxes: update 'tmp_sFlux'
        eval(store_fun);                                                    % creates/updates tmp_sFlux 

        % 90% flow routing    
        tmp_Q9_cur    = 0.9.*(flux_pn(t)-PS(tmp_sFlux(1),x1,flux_pn(t))+...
                                PERC(tmp_sFlux(1),x1)).*uh_half;            % find how the runoff of this time step will be distributed in time
        tmp_Q9_old    = tmp_Q9_old + tmp_Q9_cur;                            % update the 'still-to-flow-out' vector
        flux_q9(t)    = tmp_Q9_old(1);                                      % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_Q9_old    = circshift(tmp_Q9_old,-1);                           % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_Q9_old(end) = 0;
    
        % 10% flow routing    
        tmp_Q1_cur    = 0.1.*(flux_pn(t)-PS(tmp_sFlux(1),x1,flux_pn(t))+...
                                PERC(tmp_sFlux(1),x1)).*uh_full;            % find how the runoff of this time step will be distributed in time
        tmp_Q1_old    = tmp_Q1_old + tmp_Q1_cur;                            % update the 'still-to-flow-out' vector
        flux_q1(t)    = tmp_Q1_old(1);                                      % the first value in 'still-to-flow-out' vector leaves the model this time step
        tmp_Q1_old    = circshift(tmp_Q1_old,-1);                           % shift the 'still-to-flow-out' vector so that the next value is now at location 1
        tmp_Q1_old(end) = 0;
    
    % --- Solve store 2 (routing) -----------------------------------
        % Create temporary store ODE that need to be solved
        tmpf_S2 = @(S2)    (flux_q9(t) + ...
                            F(3.5,S2,x3,x2) - ...
                            QR(S2,x3));                                     % store 2 function
    
        % Call the numerical scheme function to create the ODE approximations
        solve_fun2 = feval(scheme,...
                        [S2old],...
                        delta_t,...
                        tmpf_S2);
                    
        % Use the specified numerical scheme to solve storages
        [tmp_sNew(2), tmp_fval] = fzero(solve_fun2,...
                                    S2old,...
                                    fzero_options);                         % 2 store solver            
         
        % --- Check if the solver has found an acceptable solution and re-run
        % if not. The re-run uses the 'lsqnonlin' solver which is slower but 
        % more robust. It runs solver.resnorm_iterations times, with different
        % starting points for the solver on each iteration ---
        tmp_resnorm = sum(tmp_fval.^2);
        if tmp_resnorm > solver.resnorm_tolerance
            [tmp_sNew(2),~,~] = rerunSolver('lsqnonlin', ...                % [tmp_sNew,tmp_resnorm,flag]
                                            lsqnonlin_options, ...          % solver options
                                            @(eq_sys) solve_fun2(...        % system of ODEs
                                                        eq_sys(1)), ...
                                            solver.resnorm_maxiter, ...     % maximum number of re-runs
                                            solver.resnorm_tolerance, ...   % convergence tolerance
                                            tmp_sNew(2), ...                % recent estimates
                                            [S2old], ...                    % storages ate previous time step
                                            store_min(2), ...               % lower bounds
                                            store_upp(2));                  % upper bounds              
        end                                       
                    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 

    % Calculate the fluxes
    flux_ps(t)   = PS(tmp_sFlux(1),x1,flux_pn(t));
    flux_es(t)   = ES(tmp_sFlux(1),x1,flux_en(t));
    flux_perc(t) = PERC(tmp_sFlux(1),x1);
    flux_qr(t)   = QR(tmp_sFlux(2),x3);
    flux_fr(t)   = F(3.5,tmp_sFlux(2),x3,x2);
    flux_fq(t)   = F(3.5,tmp_sFlux(2),x3,x2);
    flux_qt(t)   = flux_qr(t) + max(flux_q1(t) + flux_fq(t),0);
    
    % Update the stores
    store_S1(t) = S1old + (flux_ps(t) - flux_es(t) - flux_perc(t)) * delta_t;
    store_S2(t) = S2old + (flux_q9(t) + flux_fr(t) - flux_qr(t)) * delta_t;    
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the funcion and should NOT be renamed
    fluxOutput.Ea     = (flux_es + flux_ef)* delta_t;
    fluxOutput.Q      = flux_qt * delta_t;
    
    % --- Fluxes internal to the model ---
    fluxInternal.pn   = flux_pn * delta_t;
    fluxInternal.en   = flux_en * delta_t;
    fluxInternal.ef   = flux_ef * delta_t;
    fluxInternal.ps   = flux_ps * delta_t;
    fluxInternal.es   = flux_es * delta_t;
    fluxInternal.perc = flux_perc * delta_t;
    fluxInternal.q9   = flux_q9 * delta_t;
    fluxInternal.q1   = flux_q1 * delta_t;
    fluxInternal.fr   = flux_fr * delta_t;
    fluxInternal.fq   = flux_fq * delta_t;
    fluxInternal.qr   = flux_qr * delta_t;

    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;

% Check water balance
if nargout == 4

waterBalance =  sum(P).*delta_t - ...                                       % Precipitation
                sum(fluxOutput.Ea) - ...                                    % Evaporation
                sum(fluxOutput.Q) - ...                                     % Runoff
                (store_S1(end)-S10) - ...                                   % Store 1 change
                (store_S2(end)-S20) - ...                                   % Store 2 change
                (sum(tmp_Q9_old) + sum(tmp_Q1_old)).*delta_t + ...          % Still being routed
                (sum(flux_fr).*delta_t +sum(max(0,flux_q1+flux_fq).*delta_t)-sum(flux_q1).*delta_t);   % Water exchange 

    disp(['Total P  = ',num2str(sum(P).*delta_t),'mm.'])
    disp(['Total Peff  = ',num2str(sum(flux_pn)),'mm.'])
    disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
    disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
    disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
    disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
    disp(['On route = ',num2str((sum(tmp_Q9_old)+sum(tmp_Q1_old)).*delta_t),'mm.'])
    disp(['Water exchange = ',num2str(sum(flux_fr)+sum(flux_fq)),'mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a) + sum(routing)) - exchange - delta S = ',num2str(waterBalance)]);
end





