function [ theta ] = m_45_prms_18p_7s_parameter_ranges( )
%m_45_prms_18p_7s_parameter_ranges Provides parameter ranges for calibration
%   of the PRMS model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Leavesley, G. H., R. Lichty, B. Troutman, and L. Saindon (1983), 
% Precipitation-Runo Modeling System: User's Manual. U.S. Geological 
% Survey, Water-Resources Investigations Report 83-4238, 207
%
% Markstrom, S. L., S. Regan, L. E. Hay, R. J. Viger, R. M. T. Webb, R. A. 
% Payn, and J. H. LaFontaine (2015), PRMS-IV, the Precipitation-Runoff 
% Modeling System, Version 4. In U.S. Geological Survey Techniques and
% Methods, book 6, chap. B7, 158. doi: http://dx.doi.org/10.3133/tm6B7

theta = [
-3,  5;         % tt, Temperature threshold for snowfall and melt [oC]
 0, 20;         % ddf,  Degree-day factor for snowmelt [mm/oC/d]
 0,  1;         % alpha, Fraction of rainfall on soil moisture going to interception [-] 
 0,  1;         % beta, Fraction of catchment where rain goes to soil moisture [-]
 0,  5;         % stor, Maximum interception capcity [mm]
 0, 50;         % retip, Maximum impervious area storage [mm]
 0,  1;         % fscn, Fraction of SCX where SCN is located [-]
 0,  1;         % scx, Maximum contributing fraction area to saturation excess flow [-]
 0.005, 0.995;  % flz, Fraction of total soil moisture that is the lower zone [-]
 1, 2000;       % stot, Total soil moisture storage [mm]: REMX+SMAX
 0, 20;         % cgw, Constant drainage to deep groundwater [mm/d]
 1, 300;        % resmax, Maximum flow routing reservoir storage (used for scaling only, there is no overflow) [mm]
 0,  1;         % k1, Groundwater drainage coefficient [d-1]
 1,  5;         % k2, Groundwater drainage non-linearity [-]
 0,  1;         % k3, Interflow coefficient 1 [d-1]
 0,  1;         % k4, Interflow coefficient 2 [mm-1 d-1]
 0,  1;         % k5, Baseflow coefficient [d-1]
 0,  1];        % k6, Groundwater sink coefficient [d-1]
