function [ theta ] = m_18_simhyd_7p_3s_parameter_ranges( )
%m_18_simhyd_7p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the SIMHYD model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Chiew, F. H. S., Peel, M. C., & Western, A. W. (2002). Application and 
% testing of the simple rainfall-runoff model SIMHYD. In V. P. Singh & D. 
% K. Frevert (Eds.), Mathematical Models of Small Watershed Hydrology (pp. 
% 335–367). Chelsea, Michigan, USA: Water Resources Publications LLC, USA.

theta = [0,     5;      % INSC, Maximum interception capacity, [mm]
         0,     600;    % COEFF, Maximum infiltration loss parameter, [mm]
         0,     15;     % SQ, Infiltration loss exponent, [-]
         1,     2000;   % SMSC, Maximum soil moisture capacity, [mm]
         0,     1;      % SUB, Proportionality constant, [-]
         0,     1;      % CRAK, Proportionality constant, [-]
         0,     1];     % K, Slow flow time scale, [d-1]  
