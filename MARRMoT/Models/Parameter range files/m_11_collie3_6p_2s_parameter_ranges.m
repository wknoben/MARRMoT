function [ theta ] = m_11_collie3_6p_2s_parameter_ranges( )
%m_11_collie3_6p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the Collie River Basin 2 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Jothityangkoon, C., M. Sivapalan, and D. Farmer (2001), Process controls
% of water balance variability in a large semi-arid catchment: downward 
% approach to hydrological model development. Journal of Hydrology, 254,
% 174198. doi: 10.1016/S0022-1694(01)00496-6.

theta = [1   , 2000;      % Smax [mm]
         0.05, 0.95;      % fc as fraction of Smax
         0   , 1 ;        % a, subsurface runoff coefficient [d-1]
         0.05, 0.95;      % M, fraction forest cover [-]
         1   , 5          % b, flow non-linearity [-]
         0,   1];         % lambda, flow distribution [-]
