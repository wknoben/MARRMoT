function [ theta ] = m_43_gsmsocont_12p_6s_parameter_ranges( )
%m_43_gsmsocont_12p_6s_parameter_ranges Provides parameter ranges for calibration
%   of the 6-store GSM-SOCONT model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Schaefli, B., Hingray, B., Niggli, M., & Musy, a. (2005). A conceptual 
% glacio-hydrological model for high mountainous catchments. Hydrology and 
% Earth System Sciences Discussions, 2(1), 73–117. 
% http://doi.org/10.5194/hessd-2-73-2005

theta = [
            0, 1;       % fice,   Fraction of catchment covered by glacier [-]
            -3, 5;      % t0,     Threshold temperature for snowfall [oC]
            0, 20;      % asnow,  Degree-day factor for snow melt [mm/oC/d]
            -3, 3;      % tm,     Threshold temperature for snow melt [oC]
            0, 1;       % ks,     Runoff coeficient for snow melt on glacier [d-1]
            0, 20;      % aice,   Degree-day factor for ice melt [mm/oC/d]
            0, 1;       % ki,     Runoff coeficient for ice melt on glacier [d-1]
            1, 2000;    % a,      Maximum soil moisture storage [mm]
            0, 10;      % x,      Evaporation non-linearity [-]
            0, 5;       % y,      Infiltration non-linearity [-]
            0, 1;       % ksl,    Runoff coefficient for baseflow [d-1]
            0, 1];      % beta,   Runoff coefficient for quick flow [mm^(4/3)/d]
