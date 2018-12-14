function [ theta ] = m_44_echo_16p_6s_parameter_ranges( )
%m_44_echo_15p_6s_parameter_ranges Provides parameter ranges for calibration
%   of the 6-store ECHO model.
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
%
% Schaefli, B., Nicotina, L., Imfeld, C., Da Ronco, P., Bertuzzo, E., & 
% Rinaldo, A. (2014). SEHR-ECHO v1.0: A spatially explicit hydrologic 
% response model for ecohydrologic applications. Geoscientific Model 
% Development, 7(6), 2733–2746. http://doi.org/10.5194/gmd-7-2733-2014

theta = [    0, 5;      % rho, Maximum interception storage [mm]
            -3, 5;      % ts, Threshold temperature for snowfall [oC]
            -3, 3;      % tm, Threshold temperature for snowmelt [oC]
             0, 20;     % as, Degree-day factor [mm/oC/d]
             0, 1;      % af, Refreezing reduction factor [-]
             0, 2;      % gmax, Maximum melt due to ground-heat flux [mm/d]
             0, 1;      % the, Water-holding capacity of snow [-]
             0, 200;    % phi, Maximum infiltration rate [mm/d]
             1, 2000;   % smax, Maximum soil moisture storage [mm]
             0.05, 0.95;% fsm, Plant stress point as a fraction of Smax [-]
             0.05, 0.95;% fsw, Wilting point as fraction of sm [-]
             0, 1;      % ksat,  Runoff rate from soil moisture [d-1]
             0, 5;      % c, Runoff non-linearity from soil moisture [-]
             0, 20;     % lmax, Groundwater flux [mm/d]
             0, 1;      % kf, Runoff coefficient [d-1]
             0, 1];     % ks, Runoff coefficient [d-1]
             
            
