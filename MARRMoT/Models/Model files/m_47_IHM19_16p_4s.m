function [ fluxOutput, fluxInternal, storeInternal, waterBalance ] = ...
    m_47_IHM19_16p_4s( fluxInput, storeInitial, theta, solver )
% Hydrologic conceptual model: Forellenbach
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Brandes, C. (2020). Erstellung eines Konzeptionellen
% Hochwasserabfluss-modells für das Einzugsgebiet des Forellenbachs, NP
% Bayerischer Wald. MSc Thesis. Technische Univeristät Dresden, Germany.
%
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

% Data: adjust the units so that all fluxes (rates) inside this model
% function are calculated in [mm/d]
P     = fluxInput.precip./delta_t;          % [mm/delta_t] -> [mm/d]
Ep    = fluxInput.pet./delta_t;             % [mm/delta_t] -> [mm/d]
T     = fluxInput.temp;
t_end = length(P);

% Parameters
% [name in documentation] = theta(order in which specified in parameter file)

SIMAX   = theta(1);                         % maximum interception storage [mm]
A       = theta(2);                         % splitting coeffcient for excess precipitation [-]
FF      = theta(3);                         % forest fraction [-]

SMPMAX  = theta(4);                         % maximum storage macropores [mm]
CQMP    = theta(5);                         % runoff time parameter (fast/slow runnoff) first soil layer [1/d]
XQMP    = theta(6);                         % runoff scale parameter first soil layer [-]

SS1MAX  = theta(7);                         % maximum soil moisture storage first soil layer [mm]
FCCS1   = theta(8);                         % field capacity coefficient fist soil layer [-]
CFS1    = theta(9);                         % maximum infiltration rate first soil layer [-]
XFS1    = theta(10);                        % infiltration loss exponent first soil layer [-]
CQS1    = theta(11);                        % runoff time parameter for (fast/slow runnoff) first soil layer [1/d]
XQS1    = theta(12);                        % runoff scale parameter first soil layer [-]

SS2MAX  = theta(13);                        % maximum soil moisture storage second soil layer [mm]
CQS2    = theta(14);                        % runoff time parameter for (fast/slow runnoff) second soil layer [1/d]
XQS2    = theta(15);                        % runoff scale parameter second soil layer [-]

D       = theta(16);                        % Flow delay before surface runoff [d]

%INITIALISE MODEL STORES
S10     = SIMAX*0.8;                   % Initial interception storage
S20     = SMPMAX*0.8;                  % Initial macropore storage 
S30     = SS1MAX*0.8;                  % Initial soil moisture storage 1 
S40     = SS2MAX*0.8;                  % Initial soil moisture storage 2
% 
% S10     = storeInitial(1);                  % Initial interception storage
% S20     = storeInitial(2);                  % Initial macropore storage 
% S30     = storeInitial(3);                  % Initial soil moisture storage 1
% S40     = storeInitial(4);                  % Initial soil moisture storage 2

%%DEFINE STORE BOUNDARIES
store_min = [0,0,0,0];                      % lower bounds of stores
store_upp = [];                             % optional higher bounds

%%INITIALISE STORAGE VECTORS
store_S1 = zeros(1,t_end);
store_S2 = zeros(1,t_end);
store_S3 = zeros(1,t_end);
store_S4 = zeros(1,t_end);

%%INITIALISE FLUX VECTORS
flux_ei      = zeros(1,t_end);
flux_pex     = zeros(1,t_end);
flux_pexmp   = zeros(1,t_end);
flux_pexs1   = zeros(1,t_end);
flux_fmp     = zeros(1,t_end);
flux_qexmp   = zeros(1,t_end);
flux_qmp     = zeros(1,t_end);
flux_pqexs1  = zeros(1,t_end);
flux_fs1     = zeros(1,t_end);
flux_etas1   = zeros(1,t_end);
flux_qs1     = zeros(1,t_end);
flux_q0      = zeros(1,t_end);
flux_q0r     = zeros(1,t_end);
flux_qmps1   = zeros(1,t_end);
flux_pc      = zeros(1,t_end);
flux_qh      = zeros(1,t_end);
flux_qs2     = zeros(1,t_end);


qgw          = cell(1,t_end);
qgw(1,:)     = {0.0195};
flux_qgw     = cell2mat(qgw);

%%PREPARE UNIT HYDROGRAPHS
[~,uh_full] = uh_4_full(1,D,delta_t);


%%INITIALISE ROUTING VECTORS
tmp_q0r_old  = zeros(1,length(uh_full)); 


%% 1. ODEs
% Given in the documentation.

%% 2. Constitutive functions
% Given in the documentation

%% 3. Specify and smooth model functions
% Store numbering:
% S1. Interception storage SI
% S2. macropore storage SMP
% S3. soil moisture storage layer 1 SS1
% S4. soil moisture storage layer 2 SS2

% EI(SI,Ep,dt): Evaporation at the potential rate, when S>0
EI = evap_1;

% PEX(P,SI,SIMAX): Interception excess when maximum capacity is reached
PEX = interception_1;

% PEXS1(A,PEX): Split flow A
PEXS1 = split_1;

% PEXMP(A,PEX): Split flow (1-A)
PEXMP = split_2;

% FMP(SMP,SMPMAX,PEXMP): Infiltration into storage is equal the inflow when 
%                        current storage is under the maximum storage 
FMP = infiltration_8;

% QEXMP(PEXMP,FMP): Saturation excess flow from macropore storage
QEXMP = subtraction_2;

% QMP(CQMP,XQMP,SMP,dt): Non-linear interflow(variant),
%                                      
QMP = interflow_3;

% PQEXS1(PEXS1, QEXMP): Flow to Soil Moisture storage
PQEXS1 = addition_2;

% FS1(CFS1,XFS1,SS1,SS1MAX,PQEXS1): Infiltration as exponentially declining
%                                  based on relative storage
FS1 = infiltration_7;

% ETAS1(FF,FCCS1,SS1,SS1MAX,Ep,dt): evapotranspiration for Forested area
ETAS1 = evap_23;

% QS1(CQS1,FCCS1,XQS1,SS1,SS1MAX,dt): Non-linear interflow(variant),
%                                      threshhold
QS1 = interflow_12;

% Q0(PQEXS1,FS1): surface runoff
Q0 = subtraction_2;

% QMPS1(QMP, QS1): runoff from upper soil layer
QMPS1 = addition_2;

% PC(SS2,SS2MAX,QMPS1): Percolation into storage is equal the inflow when 
%                        current storage is under the maximum storage 

PC = infiltration_8;

% QH(QMPS1, PC): Interflow between soil layers
QH = subtraction_2;

% QS2(CQS2,XQS2,SS2,dt): Non-linear interflow(variant)
%                                      
QS2 = interflow_3;


%% 4. Determine numerical scheme and solver settings
% Function name of the numerical scheme
scheme            = solver.name;                                            

% Define which storage values should be used to update fluxes
[~,store_fun]     = feval(scheme,storeInitial,delta_t);                     % storeInitial is used to find the number of stores, actual values are irrelevant here

% Root-finding options
fsolve_options = optimoptions('fsolve','Display','none',...                 % [1+n stores] settings of the root finding method
                              'JacobPattern', [1,0,0,0;
                                               1,1,0,0;
                                               1,1,1,0;
                                               0,1,1,1]);                   % Specify the Jacobian pattern                                             
                                      
lsqnonlin_options = optimoptions('lsqnonlin',...                            % lsqnonlin settings for cases when fsolve fails
                                 'Display','none',...
                                 'JacobPattern', [1,0,0,0;
                                                  1,1,0,0;
                                                  1,1,1,0;
                                                  0,1,1,1],...
                                 'MaxFunEvals',1000);
%% 5. Solve the system for the full time series
for t = 1:t_end

% Model setup -------------------------------------------------------------
    % Determine the old storages
    if t == 1; S1old = S10; else; S1old = store_S1(t-1); end                % store 1 at t-1
    if t == 1; S2old = S20; else; S2old = store_S2(t-1); end                % store 2 at t-1
    if t == 1; S3old = S30; else; S3old = store_S3(t-1); end                % store 3 at t-1
    if t == 1; S4old = S40; else; S4old = store_S4(t-1); end                % store 4 at t-1
    
    % Create temporary store ODE's that need to be solved
    tmpf_S1 = ...
              @(S1,S2,S3,S4)...                                             % Change in S1 depends on ...
              (P(t)-...                                                     % Precipitation                      -                                                
               EI(S1,Ep(t),delta_t)-...                                     % Interception Evaporation from S1   -
               PEX(P(t),S1,SIMAX)...
               );                                         % Excess Precipitation                       
                    
    tmpf_S2 = ...
             @(S1,S2,S3,S4)...                                              % Change in S2 depends on ...
             (FMP(S2,SMPMAX,PEXMP(A,PEX(...
                                        P(t),S1,SIMAX...
                                        )...
                                  )...
                  )-...                                                     % Macropore Infiltration-
              QMP(CQMP,XQMP,S2,delta_t));                                   % Macropore runoff
                    
    tmpf_S3 = ...
            @(S1,S2,S3,S4)...                                               % Change in S3 depends on ...  
              (...
               FS1(...
                   CFS1,XFS1,S3,SS1MAX,...
                   PQEXS1(...
                          PEXS1(...
                                A,PEX(...
                                      P(t),S1,SIMAX...
                                      )...
                                 ),...
                          QEXMP(...
                                 PEXMP(...
                                       A,PEX(...
                                             P(t),S1,SIMAX...
                                             )...
                                       ),...
                                  FMP(...
                                      S2,SMPMAX,PEXMP(...
                                                      A,PEX(P(t),S1,SIMAX)...
                                                      )...
                                      )...
                                )...
                           )... 
                    )-...                                                   % Soil Infiltration into S3       -
               ETAS1(FF,FCCS1,S3,SS1MAX,Ep(t),delta_t)-...                  % Actual Evapotranspiration from S3  -
               QS1(CQS1,FCCS1,XQS1,S3,SS1MAX,delta_t));                     % Soil runoff from S3
                    
    tmpf_S4 = ...
              @(S1,S2,S3,S4)...                                             % Change in S4 depends on ...
               (PC(S4,SS2MAX,QMPS1(...
                                   QMP(CQMP,XQMP,S2,delta_t),...
                                   QS1(CQS1,FCCS1,XQS1,S3,SS1MAX,delta_t)...
                                   )...
                   )-...                                                    % Percolation into S4-                -
                QS2(CQS2,XQS2,S4,delta_t)...                                % Soil runoff from S4
                );                                
                    
    % Call the numerical scheme function to create the ODE approximations. 
    % This returns a new anonymous function that we solve in the next step.
    
    solve_fun = feval(scheme,...                          % numerical approximation method
                      [S1old,S2old,S3old,S4old],...       % Store values at t-1
                      delta_t,...                         % time step size
                      tmpf_S1,tmpf_S2,tmpf_S3,tmpf_S4);   % anonymous functions of ODEs    

% Model solving -----------------------------------------------------------            
    % --- Use the specified numerical scheme to solve storages ---
    [tmp_sNew,tmp_fval] = fsolve(@(eq_sys) solve_fun(...                    % call solve_fun as anonymous function
                        eq_sys(1),eq_sys(2),eq_sys(3),eq_sys(4)),...        % system of storage equations
                        [S1old,S2old,S3old,S4old],...                       % storage values on previous time step
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
                                          eq_sys(1),eq_sys(2), ...
                                          eq_sys(3),eq_sys(4)), ...
                                        solver.resnorm_maxiter, ...         % maximum number of re-runs
                                        solver.resnorm_tolerance, ...       % convergence tolerance
                                        tmp_sNew, ...                       % recent estimates
                                        [S1old,S2old,S3old,S4old], ...      % storages ate previous time step
                                        store_min, ...                      % lower bounds
                                        store_upp);                         % upper bounds              
    end
    
% Model states and fluxes -------------------------------------------------    
    % Find the storages needed to update fluxes: update 'tmp_sFlux'
    eval(store_fun);                                                        % creates/updates tmp_sFlux 
    
    % Calculate the fluxes
    flux_ei(t)      = EI(tmp_sFlux(1),Ep(t),delta_t);
    flux_pex(t)     = PEX(P(t),tmp_sFlux(1),SIMAX);
    flux_pexmp(t)   = PEXMP(A,flux_pex(t));
    flux_pexs1(t)   = PEXS1(A,flux_pex(t));
    flux_fmp(t)     = FMP(tmp_sFlux(2),SMPMAX,flux_pexmp(t));
    flux_qexmp(t)   = QEXMP(flux_pexmp(t),flux_fmp(t));
    flux_qmp(t)     = QMP(CQMP,XQMP,tmp_sFlux(2),delta_t);
    flux_pqexs1(t)  = PQEXS1(flux_pexs1(t),flux_qexmp(t));
    flux_fs1(t)     = FS1(CFS1,XFS1,tmp_sFlux(3),SS1MAX,flux_pqexs1(t));
    flux_etas1(t)   = ETAS1(FF,FCCS1,tmp_sFlux(3),SS1MAX,Ep(t),delta_t);  
    flux_qs1(t)     = QS1(CQS1,FCCS1,XQS1,tmp_sFlux(3),SS1MAX,delta_t);
    flux_q0(t)      = Q0(flux_pqexs1(t),flux_fs1(t));
    flux_qmps1(t)   = QMPS1(flux_qmp(t), flux_qs1(t));
    flux_pc(t)      = PC(tmp_sFlux(4),SS2MAX,flux_qmps1(t));
    flux_qh(t)      = QH(flux_qmps1(t), flux_pc(t));
    flux_qs2(t)     = QS2(CQS2,XQS2,tmp_sFlux(4),delta_t);
  
    
    
    % Update the stores
    store_S1(t) = S1old + (P(t) - flux_ei(t) - flux_pex(t)) * delta_t;
    store_S2(t) = S2old + (flux_fmp(t)- flux_qmp(t)) * delta_t;
    store_S3(t) = S3old + (flux_fs1(t)- flux_etas1(t) - flux_qs1(t)) * delta_t;
    store_S4(t) = S4old + (flux_pc(t) - flux_qs2(t)) * delta_t;
    
    % Routing -----------------------------------------------------------------    
    % surface runoff Q0. Apply a pre-determined (line 82)
    % triangular Unit Hydrograph routing scheme to find lagged flow Qt.
    tmp_q0r_cur      = flux_q0(t).*uh_full;     % find how the runoff of this time step will be distributed in time
    tmp_q0r_old      = tmp_q0r_old + tmp_q0r_cur;                              % update the 'still-to-flow-out' vector
    flux_q0r(t)      = tmp_q0r_old(1);                                        % the first value in 'still-to-flow-out' vector leaves the model this time step
    tmp_q0r_old      = circshift(tmp_q0r_old,-1);                             % shift the 'still-to-flow-out' vector so that the next value is now at location 1
    tmp_q0r_old(end) = 0;                                                    % remove the now last value (the one that we just routed to flux_q0r(t) and is now shifted to the end of the 'still-to-flow-out' vector)
   
    
end

%% 6. Generate outputs
    % --- Fluxes leaving the model ---
    % 'Ea' and 'Q' are used outside the function and should NOT be renamed
    fluxOutput.Ea     = (flux_ei + flux_etas1) .* delta_t;
    fluxOutput.Q      = (flux_q0r + flux_qh + flux_qs2) .* delta_t + flux_qgw;
    
    % --- Fluxes internal to the model ---
    fluxInternal.ei     = flux_ei       * delta_t;
    fluxInternal.pex    = flux_pex      * delta_t;
    fluxInternal.pexmp  = flux_pexmp    * delta_t;
    fluxInternal.pexs1  = flux_pexs1    * delta_t;
    fluxInternal.fmp    = flux_fmp      * delta_t;
    fluxInternal.qexmp  = flux_qexmp    * delta_t;
    fluxInternal.qmp    = flux_qmp      * delta_t;
    fluxInternal.pqexs1 = flux_pqexs1   * delta_t;
    fluxInternal.fs1    = flux_fs1      * delta_t;
    fluxInternal.etas1  = flux_etas1    * delta_t;
    fluxInternal.qs1    = flux_qs1      * delta_t;
    fluxInternal.q0     = flux_q0       * delta_t;
    fluxInternal.q0r    = flux_q0r      * delta_t;
    fluxInternal.qmps1  = flux_qmps1    * delta_t;
    fluxInternal.pc     = flux_pc       * delta_t;
    fluxInternal.qh     = flux_qh       * delta_t;
    fluxInternal.qs2    = flux_qs2      * delta_t;
    fluxInternal.qgw    = flux_qgw;       
   
       
    % --- Stores ---
    storeInternal.S1  = store_S1;
    storeInternal.S2  = store_S2;
    storeInternal.S3  = store_S3;
    storeInternal.S4  = store_S4;
    
% % Check water balance
if nargout == 4

    % Code for manual water balance calculation
    waterBalance = sum(P.*delta_t) + ...                % Precipitation
                   sum(fluxInternal.qgw) -...           % Groundwater runoff
                   sum(fluxOutput.Ea) - ...             % Evaporation
                   sum(fluxOutput.Q) - ...              % Outflow
                   (store_S1(end)-S10) - ...            % Storage change
                   (store_S2(end)-S20)-...
                   (store_S3(end)-S30)-...
                   (store_S4(end)-S40) -...
                    sum(tmp_q0r_old.*delta_t);                     % Still in routing vector; 
    
     disp(['Total P  = ',num2str(sum(P.*delta_t)),'mm.'])
     disp(['Total Q  = ',num2str(sum(fluxOutput.Q)),'mm.'])
     disp(['Total E  = ',num2str(sum(fluxOutput.Ea)),'mm.'])
     disp(['Delta S1 = ',num2str((store_S1(end)-S10)),'mm.'])
     disp(['Delta S2 = ',num2str((store_S2(end)-S20)),'mm.'])
     disp(['Delta S3 = ',num2str((store_S3(end)-S30)),'mm.'])
     disp(['Delta S4 = ',num2str((store_S4(end)-S40)),'mm.'])
     disp(['On route = ',num2str(sum(tmp_q0r_old.*delta_t)),'mm.'])
     disp(['Water balance = sum(P) + sum(GW) - sum(Q) - sum(E_a) - sum(routing) - delta S = ',num2str(waterBalance)]);  
end




 




