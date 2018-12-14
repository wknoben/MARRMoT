function [ theta ] = m_39_mcrm_16p_5s_parameter_ranges( )
%m_39_mcrm_16p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the 5-store Midlands Catchment Runoff model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Moore, R. J., & Bell, V. A. (2001). Comparison of rainfall-runoff models 
% for flood forecasting. Part 1: Literature review of models. Bristol: 
% Environment Agency.

theta = [0   , 5;       % smax, Maximum interception storage [mm]
         0.01, 0.99;    % cmax, Maximum fraction of area contributing to rapid runoff [-]
         0.01, 0.99;    % ct, Fraction of cmax that is the minimum contributing area [-]
         0   , 2;       % c1, Shape parameter for rapid flow distribution [mm-1]
         0   , 1;       % ce, Shape parameter for evaporation [mm-1]
         1   , 2000;    % dsurp, Threshold for direct runoff [mm]
         0   , 1;       % kd, Direct runoff time parameter [d-1]
         1   , 5;       % gamd, Direct runoff flow non-linearity [-]
         0   , 20;      % qpmax, Maximum percolation rate [mm/d]
         0   , 1;       % kg, Groundwater time parameter [d-1]
         1   , 120;     % tau, Routing delay [d]
         1   , 300;     % sbf, Maximum routing store depth [mm]
         0   , 1;       % kcr, Channel flow time parameter [d-1]
         1   , 5;       % gamcr, Channel flow non-linearity [-]
         0   , 1;       % kor, Out-of-bank flow time parameter [d-1]
         1   , 5];      % gamor, Out-of-bank flow non-linearity [-]         