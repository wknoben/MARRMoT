function [ theta ] = m_15_plateau_8p_2s_parameter_ranges( )
%m_15_plateau_8p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store plateau model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Savenije, H. H. G. (2010). “Topography driven conceptual modelling 
% (FLEX-Topo).” Hydrology and Earth System Sciences, 14(12), 2681–2692. 
% https://doi.org/10.5194/hess-14-2681-2010

theta = [0   , 200;    % Fmax, maximum infiltration rate [mm/d]
         0   , 5;      % Dp, interception capacity [mm]
         1   , 2000;   % SUmax, soil misture depth [mm]
         0.05,    0.95;% Swp, wilting point as fraction of Sumax [-]
         0   , 1;      % p, coefficient for moisture constrained evaporation [-]
         1   , 120;    % tp, time delay for routing [d]
         0   , 4;      % c, capillary rise [mm/d]
         0   , 1];     % kp, base flow time parameter [d-1]
