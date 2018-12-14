function [ theta ] = m_26_flexi_10p_4s_parameter_ranges( )
%m_26_flexi_10p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store Flex-I model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Fenicia, F., McDonnell, J. J., & Savenije, H. H. G. (2008). Learning from
% model improvement: On the contribution of complementary data to process 
% understanding. Water Resources Research, 44(6), 1–13. 
% http://doi.org/10.1029/2007WR006386

theta = [1,2000;        % Maximum soil moisture storage [mm]
         0, 10;         % Unsaturated zone shape parameter [-]
         0, 1;          % Fast/slow runoff distribution parameter [-]
         0, 20;         % Maximum percolation rate [mm/d]
         0.05, 0.95;    % Wilting point as fraction of s1max [-]
         1, 5;          % Flow delay before fast runoff [d]
         1, 15;         % Flow delay before slow runoff [d]
         0, 1;          % Fast runoff coefficient [d-1]
         0, 1;          % Slow runoff coefficient [d-1]
         0, 5];         % Maximum interception storage [mm]