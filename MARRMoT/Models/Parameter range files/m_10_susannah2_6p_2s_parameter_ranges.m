function [ theta ] = m_10_susannah2_6p_2s_parameter_ranges( )
%m_10_susannah2_6p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store Susannah Brook v2 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Son, K., & Sivapalan, M. (2007). Improving model structure and reducing 
% parameter uncertainty in conceptual water balance models through the use 
% of auxiliary data. Water Resources Research, 43(1). 
% https://doi.org/10.1029/2006WR005032

theta = [1   , 2000;    % Sb, Maximum soil moisture storage [mm]
         0.05, 0.95;    % phi, Porosity [-]
         0.01, 1.00;    % Sfc,Wiliting point as fraction of sb [-]
         0   , 1 ;      % r, Fraction of runoff coefficient [-]
         0   , 1 ;      % c, Runoff coefficient [d-1] (should be > 0)
         1   , 5];      % d, Runoff coefficient [-] (should be > 0)
