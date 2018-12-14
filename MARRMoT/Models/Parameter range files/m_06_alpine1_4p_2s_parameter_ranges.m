function [ theta ] = m_06_alpine1_4p_2s_parameter_ranges( )
%m_06_alpine1_4p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store Alpine model v1.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Eder, G., Sivapalan, M., & Nachtnebel, H. P. (2003). Modelling water 
% balances in an Alpine catchment through exploitation of emergent 
% properties over changing time scales. Hydrological Processes, 17(11), 
% 2125–2149. http://doi.org/10.1002/hyp.1325

theta = [-3  , 5        % tt [celsius]
         0   , 20;      % degree-day-factor
         1   , 2000;    % Smax [mm]
         0   , 1];      % time delay of subsurface flow [d-1]
