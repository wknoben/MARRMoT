function [ theta ] = m_47_IHM19_16p_4s_parameter_ranges( )
%m_47_IHM19_16p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store IHM19 model.
%
% Copyright (C) 2021 Clara Brandes
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Brandes, C. (2020). Erstellung eines Konzeptionellen
% Hochwasserabfluss-modells für das Einzugsgebiet des Forellenbachs, NP
% Bayerischer Wald. MSc Thesis. Technische Univeristät Dresden, Germany.

theta = [   2,     5;      % 01 SIMAX, maximum interception storage [mm]
            0.9,   1;      % 02 A, splitting coeffcient for excess precipitation [-]
            0.4,   0.95;   % 03 FF, forest fraction [-]
            0.05,  5;      % 04 SMPMAX, maximum storage macropores [mm]
            0,     1;      % 05 CQMP, runoff time parameter (fast/slow runnoff) first soil layer [1/d]
            1,     5;      % 06 XQMP, runoff scale parameter first soil layer [-]     
          400,   600;      % 07 SS1MAX, maximum soil moisture storage first soil layer [mm]
            0.3,   0.7;    % 08 FCCS1, field capacity as fraction of maximum storage first soil layer [-]
            0,  1000;      % 09 CFS1, maximum infiltration rate first soil layer [mm/d]
            0,    15;      % 10 XFS1, infiltration loss exponent first soil layer [-]
            0,     1;      % 11 CQS1, runoff time parameter for (fast/slow runnoff) first soil layer [1/d]
            1,     5;      % 12 XQS1, runoff scale parameter first soil layer [-]
          300,   500;      % 13 SS2MAX, maximum soil moisture storage second soil layer [mm]
            0,     1;      % 14 CQS2, runoff time parameter for (fast/slow runnoff) second soil layer [1/d]
            1,     5;      % 15 XQS2, runoff scale parameter second soil layer [-]
            0.01,  5];     % 16 D, Flow delay before surface runoff [d]
  
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'SIMAX, maximum interception storage [mm]' ...
                           'A, splitting coeffcient for excess precipitation [-]' ...
                           'FF, forest fraction' ...
                           'SMPMAX, maximum storage macropores [mm]' ...
                           'CQMP, runoff time parameter (fast/slow runnoff) first soil layer [1/d]' ...
                           'XQMP, runoff scale parameter first soil layer [-]'...
                           'SS1MAX, maximum soil moisture storage first soil layer [mm]'...
                           'FCCS1, field capacity coefficient fist soil layer [-]'...
                           'CFS1, maximum infiltration rate first soil layer [-]'...
                           'XFS1, infiltration loss exponent first soil layer [-]'...
                           'CQS1, runoff time parameter for (fast/slow runnoff) first soil layer [1/d]'...
                           'XQS1, runoff scale parameter first soil layer [-]'...
                           'SS2MAX, maximum soil moisture storage second soil layer [mm]'...
                           'CQS2, runoff time parameter for (fast/slow runnoff) second soil layer [1/d]'...
                           'XQS2, runoff scale parameter second soil layer [-]'...
                           'D, Flow delay before surface runoff [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp('Please note that the IHM19 parameter ranges are set for a specific experiment (see model reference) and not necessarily consistent with the parameter ranges of MARRMoT models 1-46.')
disp(txt)