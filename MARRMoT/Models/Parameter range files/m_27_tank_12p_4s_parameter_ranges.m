function [ theta ] = m_27_tank_12p_4s_parameter_ranges( )
%m_27_tank_12p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store Tank Model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Sugawara, M. (1995). Tank model. In V. P. Singh (Ed.), Computer models of
% watershed hydrology (pp. 165–214). Water Resources Publications, USA.

theta = [   0, 1;       % Time parameter for drainage 1>2 [d-1]
            0, 1;       % Time parameter for drainage 2>3 [d-1]
            0, 1;       % Time parameter for drainage 3>4 [d-1]
            0, 1;       % Time parameter for surface runoff 1 [d-1]
            0, 1;       % Fraction of a1 that is a2 [-]
            0, 1;       % Fraction of a2 that is b1 [-]   
            0, 1;       % Fraction of b1 that is c1 [-]
            0, 1;       % Fraction of c1 that is d1 [-]
            1, 2000;    % Maximum soil depth (sum of runoff thresholds) [mm]
            0.01, 0.99; % Fraction of st that consitutes threshold t2 [-]
            0.01, 0.99; % Fraction of st-t2 that is added to t2 to find threshold 1 [-] (ensures t1 > t2)
            0.01, 0.99];% Fraction of st-t1-t2 that consitutes threshold 3 [-]