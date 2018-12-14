function [ theta ] = m_04_newzealand1_6p_1s_parameter_ranges( )
%m_04_newzealand1_6p_1s_parameter_ranges Provides parameter ranges for 
% calibration of the 1-store New Zealand model v1.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Atkinson, S. E., Woods, R. A., & Sivapalan, M. (2002). Climate and
% landscape controls on water balance model complexity over changing 
% timescales. Water Resources Research, 38(12), 17–50. 
% http://doi.org/10.1029/2002WR001487

theta = [1   , 2000;   % Smax, Maximum soil moisture storage [mm]
         0.01, 1;      % sfc, Field capacity as fraction of maximum soil moisture [-]
         0.05, 0.95 ;  % m, Fraction forest [-]
         0   , 1;      % a, Subsurface runoff coefficient [d-1]
         1   , 5;      % b, Runoff non-linearity [-]
         0   , 1];     % tcbf, Baseflow runoff coefficient [d-1]
     
                        
