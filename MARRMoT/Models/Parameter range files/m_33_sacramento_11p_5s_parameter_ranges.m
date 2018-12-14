function [ theta ] = m_33_sacramento_11p_5s_parameter_ranges( )
%m_33_sacramento_11p_5s_parameter_ranges Provides parameter ranges for 
% calibration of the SACRAMENTO model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% National Weather Service (2005), II.3-SAC-SMA: Conceptualization of the 
% Sacramento Soil Moisture Accounting model.In National Weather Service 
% River Forecast System (NWSRFS) User Manual, 113
%
% Koren, V. I., Smith, M., Wang, D., & Zhang, Z. (2000). Use of soil 
% property data in the derivation of conceptual rainfall-runoff model 
% parameters. Proceedings of the 15th Conference on Hydrology, AMS, Long 
% Beach, CA, (1), 103–106.

theta = [   
0       , 1;        % pctim, Fraction impervious area [-]
1       , 2000;     % smax, Maximum total storage depth [mm]
0.005   , 0.995;    % f1, fraction of smax that is Maximum upper zone tension water storage [mm]
0.005   , 0.995;    % f2, fraction of smax-S1max that is Maximum upper zone free water storage [mm]
0       , 1;        % kuz, Interflow runoff coefficient [d-1]
0       , 7;        % rexp, Base percolation rate non-linearity factor [-]
0.005   , 0.995;    % f3, fraction of smax-S1max-S2max that is  Maximum lower zone tension water storage [mm]
0.005   , 0.995;    % f4, fraction of smax-S1max-S2max-S3max that is  Maximum lower zone primary free water storage [mm]
0       , 1;        % pfree, Fraction of percolation directed to free water stores [-]
0       , 1;        % klzp, Primary baseflow runoff coefficient [d-1]
0       , 1];       % klzs, Supplemental baseflow runoff coefficient [d-1]
