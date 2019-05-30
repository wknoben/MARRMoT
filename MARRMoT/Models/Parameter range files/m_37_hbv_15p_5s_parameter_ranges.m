function [ theta ] = m_37_hbv_15p_5s_parameter_ranges( )
%m_37_hbv_15p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the HBV-96 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Lindström, G., Johansson, B., Persson, M., Gardelin, M., & Bergström, S. 
% (1997). Development and test of the distributed HBV-96 hydrological model. 
% Journal of Hydrology, 201, 272–288. 
% https://doi.org/https://doi.org/10.1016/S0022-1694(97)00041-3

theta = [-3 ,   5;      % TT, threshold temperature for snowfall [oC] 
          0 ,   17;     % TTI, interval length of rain-snow spectrum [oC]
         -3 ,   3;      % TTM, threshold temperature for snowmelt [oC]
          0 ,   1;      % CFR, coefficient of refreezing of melted snow [-]
          0,   20;      % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
          0 ,   1;      % WHC, maximum water holding content of snow pack [-]
          0 ,   4;      % CFLUX, maximum rate of capillary rise [mm/d]
          1 ,   2000;   % FC, maximum soil moisture storage [mm]
          0.05,    0.95;% LP, wilting point as fraction of FC [-]
          0 ,   10;     % BETA, non-linearity coefficient of upper zone recharge [-]
          0 ,   1;      % K0, runoff coefficient from upper zone [d-1] 
          0 ,   4;      % ALPHA, non-linearity coefficient of runoff from upper zone [-]
          0 ,   20;     % PERC, maximum rate of percolation to lower zone [mm/d]
          0 ,   1;      % K1, runoff coefficient from lower zone [d-1]
          1 ,   120];   % MAXBAS, flow routing delay [d]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'TT, threshold temperature for snowfall [oC]' ...
                           'TTI, interval length of rain-snow spectrum [oC]' ...
                           'TTM, threshold temperature for snowmelt [oC]' ...
                           'CFR, coefficient of refreezing of melted snow [-]' ...
                           'CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]' ...
                           'WHC, maximum water holding content of snow pack [-]' ...
                           'CFLUX, maximum rate of capillary rise [mm/d]' ...
                           'FC, maximum soil moisture storage [mm]' ...
                           'LP, wilting point as fraction of FC [-]' ...
                           'BETA, non-linearity coefficient of upper zone recharge [-]' ...
                           'K0, runoff coefficient from upper zone [d-1]' ...
                           'ALPHA, non-linearity coefficient of runoff from upper zone [-]' ...
                           'PERC, maximum rate of percolation to lower zone [mm/d]' ...
                           'K1, runoff coefficient from lower zone [d-1]' ...
                           'MAXBAS, flow routing delay [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)      
