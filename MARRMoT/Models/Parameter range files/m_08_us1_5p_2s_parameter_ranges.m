function [ theta ] = m_08_us1_5p_2s_parameter_ranges( )
%m_08_us1_5p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store test United States model v1.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Bai, Y., Wagener, T., & Reed, P. (2009). A top-down framework for 
% watershed model evaluation and selection under uncertainty. Environmental
% Modelling & Software, 24(8), 901–916. 
% http://doi.org/10.1016/j.envsoft.2008.12.012

theta = [0   , 1;       % Alpha_ei, Fraction of intercepted rainfall [-]
         0.05, 0.95;    % M, Fraction forest [-]
         1   , 2000;    % Smax, Maximum soil moisture [mm]
         0.01, 1;       % fc, Field capacity as fraction of Smax [-]
         0   , 1];      % Alpha_ss, Subsurface routing delay [d-1]
