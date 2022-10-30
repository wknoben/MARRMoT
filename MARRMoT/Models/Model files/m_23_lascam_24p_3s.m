classdef m_23_lascam_24p_3s < MARRMoT_model
% Class for hydrologic conceptual model: Large-scale catchment water and 
% salt balance model element 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Model reference
% Sivapalan, M., Ruprecht, J. K., & Viney, N. R. (1996). Water and salt 
% balance modelling to predict the effects of land-use changes in forested 
% catchments. 1. Small catchment water balance model. Hydrological 
% Processes, 10(3).

    properties
        % model-specific attributes
        aux_theta     % auxiliary parameters
    end
    methods
        
        % creator method
        function obj = m_23_lascam_24p_3s()
            obj.numStores = 3;                                             % number of model stores
            obj.numFluxes = 16;                                            % number of model fluxes
            obj.numParams = 24; 

            obj.JacobPattern  = [1,1,1;
                                 1,1,1;
                                 1,1,1];                                   % Jacobian matrix of model store ODEs
                             
            obj.parRanges = [   0, 200;         % af, Catchment-scale infiltration parameter [mm/d]
                                0, 5;           % bf, Catchment-scale infiltration non-linearity parameter [-]
                                1, 2000;        % stot, Total catchment storage [mm]
                                0.01, 0.99;     % xa, Fraction of Stot that is Amax [-]
                                0.01, 0.99;     % xf, Fraction of Stot-Amx that is depth Fmax [-]
                                0.01, 0.99;     % na, Fraction of Amax that is Amin [-]
                                0, 5;           % ac, Variable contributing area scaling [-]
                                0, 10;          % bc, Variable contributing area non-linearity [-]
                                0, 5;           % ass, Subsurface saturation area scaling [-]
                                0, 10;          % bss, Subsurface saturation area non-linearity [-]
                                0, 200;         % c, Maximum infiltration rate [mm/d]
                                0, 5;           % ag, Interception base parameter [mm/d]
                                0, 1;           % bg, Interception fraction parameter [-]
                                0, 1;           % gf, F-store evaporation scaling [-]
                                0, 10;          % df, F-store evaporation non-linearity [-]
                                0, 1;           % rd, Recharge time parameter [d-1]
                                0, 1;           % ab, Groundwater flow scaling [-]
                                0.01, 200;      % bb, Groundwater flow base rate [mm/d]
                                0, 1;           % ga, A-store evaporation scaling [-]
                                0, 10;          % da, A-store evaporation non-linearity [-]
                                0.01, 200;      % aa, Subsurface storm flow rate [mm/d]
                                1, 5;           % ba, Subsurface storm flow non-linearity [-]
                                0, 1;           % gb, B-store evaporation scaling [-]
                                0, 10];         % db, B-store evaporation non-linearity [-]
            
            obj.StoreNames = {"S1", "S2" "S3"};                             % Names for the stores
            obj.FluxNames  = {"ei",   "pg",   "qse", "qie", "pc",...
                              "qsse", "qsie", "fa",  "ef",  "rf",...
                              "ea1",  "ea2",  "qa",  "ra",  "qb", "eb"};   % Names for the fluxes
            
            obj.FluxGroups.Ea = [1 9 11 12 16];                            % Index or indices of fluxes to add to Actual ET
            obj.FluxGroups.Q  = [3 4 13];                                  % Index or indices of fluxes to add to Streamflow
            
        end
        
        % INITialisation function
        function obj = init(obj)
            % parameters
            theta = obj.theta;
            stot = theta(3);     % Total catchment storage [mm]
            xa   = theta(4);     % Fraction of Stot that is Amax [-]
            xf   = theta(5);     % Fraction of Stot-Amx that is depth Fmax [-]
            na   = theta(6);     % Fraction of Amax that is Amin [-];
            
            % auxiliary parameters
            amax = xa*stot;             % Maximum contributing area depth [mm]
            fmax = xf*(stot-amax);      % Infiltration depth scaling [mm]
            bmax = (1-xf)*(stot-amax);  % Groundwater depth scaling [mm]
            amin = na*amax;             % Minimum contributing area depth [mm]
            obj.aux_theta = [amax fmax bmax amin];
            
        end
        
        % MODEL_FUN are the model governing equations in state-space formulation
        function [dS, fluxes] = model_fun(obj, S)
            % parameters
            theta = obj.theta;
            af   = theta(1);     % Catchment-scale infiltration parameter [mm/d]
            bf   = theta(2);     % Catchment-scale infiltration non-linearity parameter [-]
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
            
            % auxiliary parameters
            aux_theta = obj.aux_theta;
            amax = aux_theta(1);  % Maximum contributing area depth [mm]
            fmax = aux_theta(2);  % Infiltration depth scaling [mm]
            bmax = aux_theta(3);  % Groundwater depth scaling [mm]
            amin = aux_theta(4);  % Minimum contributing area depth [mm]
            
            % delta_t
            delta_t = obj.delta_t;
            
            % stores
            S1 = S(1);
            S2 = S(2);
            S3 = S(3);
            
            % climate input
            t = obj.t;                             % this time step
            climate_in = obj.input_climate(t,:);   % climate at this step
            P  = climate_in(1);
            Ep = climate_in(2);
            T  = climate_in(3);
            
            % fluxes functions
            tmp_phiss= area_1(ass,bss,S2,amin,amax);
            tmp_phic = area_1(ac,bc,S2,amin,amax);
            tmp_fss  = infiltration_5(af,bf,S3,bmax,S1,fmax);

            flux_pg  = interception_5(bg,ag,P);
            flux_ei  = effective_1(P,flux_pg);
            flux_qse = saturation_11(ac,bc,S2,amin,amax,flux_pg);
            flux_pc  = infiltration_4(flux_pg-flux_qse,c);
            flux_qie = effective_1(flux_pg-flux_qse,flux_pc);
            flux_qsse= saturation_12(tmp_phiss,tmp_phic,flux_pc);
            flux_fa  = infiltration_4(max(0,flux_pc*min(1,(1-tmp_phiss)/(1-tmp_phic))),tmp_fss);
            flux_qsie= effective_1(flux_pc,flux_fa+flux_qsse);
            flux_ef  = evap_19(gf,df,S1,fmax,Ep,delta_t);
            flux_rf  = recharge_3(td,S1);
            flux_ea1 = evap_1(S2,tmp_phic*Ep,delta_t) ;
            flux_ea2 = evap_19(ga,da,S2,amax,Ep,delta_t);
            flux_qa  = saturation_11(aa,ba,S2,amin,amax,1);
            flux_ra  = recharge_4(tmp_phic,tmp_fss,delta_t);
            flux_qb  = baseflow_8(bb,ab,S3,bmax);
            flux_eb  = evap_19(gb,db,S3,bmax,Ep,delta_t);

            % stores ODEs
            dS1 = flux_fa - flux_ef - flux_rf;
            dS2 = flux_qsse + flux_qsie + flux_qb - flux_ea1 - ...
                  flux_ea2  - flux_ra   - flux_qa;
            dS3 = flux_rf + flux_ra - flux_eb - flux_qb; 
            
            % outputs
            dS = [dS1 dS2 dS3];
            fluxes = [flux_ei, flux_pg,   flux_qse,  flux_qie,...
                      flux_pc, flux_qsse, flux_qsie, flux_fa,...
                      flux_ef, flux_rf,   flux_ea1,  flux_ea2,...
                      flux_qa, flux_ra,   flux_qb,   flux_eb];
        end
        
        % STEP runs at the end of every timestep
        function obj = step(obj)
        end
    end
end