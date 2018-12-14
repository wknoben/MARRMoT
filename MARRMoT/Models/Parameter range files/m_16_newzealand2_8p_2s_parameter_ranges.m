function [ theta ] = m_16_newzealand2_8p_2s_parameter_ranges( )
%m_16_newzealand2_8p_2s_parameter_ranges Provides parameter ranges for 
% calibration of the 1-store New Zealand model v1.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Atkinson, S. E., Sivapalan, M., Woods, R. A., & Viney, N. R. (2003). 
% Dominant physical controls on hourly flow predictions and the role of 
% spatial variability: Mahurangi catchment, New Zealand. Advances in Water 
% Resources, 26(3), 219–235. http://doi.org/10.1016/S0309-1708(02)00183-5

theta = [0   , 5;      % Maximum interception storage [mm] 
         1   , 2000;   % Smax, Maximum soil moisture storage [mm]
         0.01, 1;      % sfc, Field capacity as fraction of maximum soil moisture [-]
         0.05, 0.95 ;  % m, Fraction forest [-]
         0   , 1;      % a, Subsurface runoff coefficient [d-1]
         1   , 5;      % b, Runoff non-linearity [-]
         0   , 1;      % tcbf, Baseflow runoff coefficient [d-1]
         1   , 120];   % Routing time delay [d]
     
                        
