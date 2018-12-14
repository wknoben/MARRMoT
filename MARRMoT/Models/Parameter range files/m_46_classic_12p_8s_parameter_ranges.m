function [ theta ] = m_46_classic_12p_8s_parameter_ranges( )
%m_46_classic_12p_8s_parameter_ranges Provides parameter ranges for calibration
%   of the 8-store CLASSIC model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Crooks, S. M., & Naden, P. S. (2007). CLASSIC: a semi-distributed 
% rainfall-runoff modelling system. Hydrology and Earth System Sciences, 
% 11(1), 516–531. http://doi.org/10.5194/hess-11-516-2007

theta = [   0, 1;       % fap, Fraction of catchment area that has permeable soils [-]
            0.01, 0.99; % fdp, Fraction of depth of permeable soil that is store Px [-]
            1, 2000;    % Depth of permeable soil [mm]
            0, 1;       % Runoff coefficient for permeable soil [d-1]
            0, 1;       % Fraction of Ps that infiltrates into semi-permeable soil [-]
            0, 1;       % Fraction of (1-fap) that is fas [-]
            0.01, 0.99; % fds, Fraction of depth of semi-permeable soil that is store Sx [-]
            1, 2000;    % Depth of semi-permeable soil [mm]
            0, 1;       % Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]
            0, 1;       % Quick runoff coefficient for semi-permeable soil [d-1]
            0, 1;       % Slow runoff coefficient for semi-permeable soil [d-1]
            0, 1];      % Runoff coefficient for impermeable soil [d-1]