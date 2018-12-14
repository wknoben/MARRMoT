function [ theta ] = m_09_susannah1_6p_2s_parameter_ranges( )
%m_09_susannah1_6p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store Susannah Brook v1-5 model.
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
         0.05,    0.95; % Sfc,Wiliting point as fraction of sb [-]
         0   , 1 ;      % M, Fraction forest [-]
         1   , 50;      % a, Runoff coefficient [d] (should be > 0)
         0.2 , 1;       % b, Runoff coefficient [-] (should be > 0)
         0   , 1];      % r, Runoff coefficient [d-1] 
