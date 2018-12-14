function [ theta ] = m_nn_example_7p_3s_parameter_ranges( )
%m_nn_example_7p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the 3-store example model, created by W. Knoben in 09-2018.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

theta = [ 0,    4;     % crate, Maximum capillary rise rate [mm/d]
          1, 2000;     % uzamx, Maximum upper zone storage [mm]
          0,   20;     % prate, Maximum percolation rate [mm/d]
          0,    1;     % klz, Lower zone runoff coefficient [d-1]
          0,    1;     % alpha, Fraction of lower zone runoff to groundwater [-]
          0,    1;     % kg, Groundwater runoff coefficient [d-1]
          1,  120];    % d, Routing delay [d]
