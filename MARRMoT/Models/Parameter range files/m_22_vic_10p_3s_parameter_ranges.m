function [ theta ] = m_22_vic_10p_3s_parameter_ranges( )
%m_22_vic_10p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the 3-store VIC model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. A., Vrugt, J. A., 
% Gupta, H. V., … Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between 
% hydrological models. Water Resources Research, 44(12), W00B02. 
% http://doi.org/10.1029/2007WR006735
%
% Liang, X., Lettenmaier, D. P., Wood, E. F., & Burges, S. J. (1994). A 
% simple hydrologically based model of land surface water and energy fluxes
% for general circulation models. Journal of Geophysical Research, 99, 
% 14415–14428.

% NOTE: the chosen range for interception capacity is [0.1, 5]. For other
% models this is [0, 5]. However, most other models assume that evaporation
% from interception occurs at the potential rate, whereas in VIC this
% evaporation is assumed to occur at a fraction S/Smax of PET. If Smax = 0,
% this leads to numerical issues. Therefor the maximum interception range
% in this model is slightly different from that in other models.


theta = [   0.1 , 5;        % ibar,     Mean interception capacity [mm]
            0   , 1;        % idelta,   Seasonal interception change as fraction of mean [-]
            1   , 365;      % ishift,   Maximum interception peak timing [-]
            1   , 2000;     % stot,    Maximum soil moisture capacity [mm]
            0.01, 0.99;     % fsm, Fraction of stot that constitutes maximum soil moisture smmax [-]
            0   ,10;        % b,        Infiltration excess shape parameter [-]
            0   , 1;        % k1,       Percolation time parameter [d-1]
            0   ,10;        % c1,       Percolation non-linearity parameter [-]
            0   ,1;         % k2,       Baseflow time parameter [d-1]
            1   , 5];       % c2,       Baseflow non-linearity parameter