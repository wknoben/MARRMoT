function [ theta ] = m_35_mopex5_12p_5s_parameter_ranges( )
%m_35_mopex5_12p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the 5-store MOPEX-5 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Ye, S., Yaeger, M., Coopersmith, E., Cheng, L., & Sivapalan, M. (2012). 
% Exploring the physical controls of regional patterns of flow duration 
% curves - Part 2: Role of seasonality, the regime curve, and associated 
% process controls. Hydrology and Earth System Sciences, 16(11), 4447–4465.
% http://doi.org/10.5194/hess-16-4447-2012

theta = [-3  , 3;       % tcrit, Snowfall & snowmelt temperature [oC]
         0   , 20;      % ddf, Degree-day factor for snowmelt [mm/oC/d]
         1   , 2000;    % Sb1, Maximum soil moisture storage [mm]
         0   , 1 ;      % tw, Groundwater leakage time [d-1]
         0   , 1 ;      % I_alpha, Intercepted fraction of Pr [-]
         1   , 365;     % I_s, Maximum Leaf Area Index timing [d]
         -10 , 0;       % tmin, Growing Season Index minimum temperature
         1   , 20;      % trange, Growing Season Index temperature range
         0   , 1 ;      % tu, Slow flow routing response time [d-1]
         0.05, 0.95;    % se, Root zone storage capacity as fraction of Sb2 [-]
         1   , 2000;    % Sb2, Root zone storage capacity [mm]
         0   , 1];      % tc, Mean residence time [d-1]
