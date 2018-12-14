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

theta = [-3 ,   5;      % TT, threshold temperature for snowfall [oC]. 
          0 ,   17;     % TTI, interval length of rain-snow spectrum [oC]
         -3 ,   3;      % TTM, threshold temperature for snowmelt [oC]
          0 ,   1;      % CFR, coefficient of refreezing of melted snow [-]
          0,   20;      % CFMAX, degree-day factor of snowmelt and refreezing [mm/oC/d]
          0 ,   1;      % WHC, maximum water holding content of snow pack [-]
          0 ,   4;      % CFLUX, maximum rate of capillary rise [mm/d]
          1 ,   2000;   % FC, maximum soil moisture storage [mm]
          0.05,    0.95;% LP, wilting point as fraction of FC [-]
          0 ,   10;     % BETA, non-linearity coefficient of upper zone recharge [-]
          0 ,   1;      % K0, runoff coefficient from upper zone [d-1], 
          0 ,   4;      % ALPHA, non-linearity coefficient of runoff from upper zone [-]
          0 ,   20;     % PERC, maximum rate of percolation to lower zone [mm/d]
          0 ,   1;      % K1, runoff coefficient from lower zone [d-1]
          1 ,   120];   % MAXBAS, flow routing delay [d]
      
      
% theta = [-3   ,   3;         x % TT [oC].          Def = 
%           0   ,   7 ;        x % TTI [oC].         Def = 2.00
%          -x   ,   x ;        x % TTM [oC].         Def = 
%           0   ,   1;         x % CFR [-].          Def = 0.05
%           0   ,   20;        x % CFMAX [mm/oC/d].  Def = 3.50
%           0   ,   0.5;       x % WHC [-].          Def = 0.10
%           0   ,   4;         x % CFLUX [mm/d].     Def = 1.00
%           1   ,   2000;      x % FC [mm].          Def = 
%           0.1 ,   1;          % LP [-].           Def = 
%           1   ,   6;          % BETA [-].         Def = 2.00
%           0.0005 ,0.3;       x % K0 [d-1].         Def = 0.0005
%           0   ,   3;          % ALPHA [-].        Def = 
%           0   ,   6;         x % PERC [mm/d].      Def = 
%           0.0005 ,0.3;       x % K1 [d-1].         Def = 
%           0   ,   x];        x % MAXBAS [d].       Def = 