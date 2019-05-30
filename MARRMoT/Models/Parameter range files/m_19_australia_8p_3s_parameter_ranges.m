function [ theta ] = m_19_australia_8p_3s_parameter_ranges( )
%m_19_australia_8p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the 3-store Australia model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Farmer, D., Sivapalan, M., & Jothityangkoon, C. (2003). Climate, soil, 
% and vegetation controls upon the variability of water balance in 
% temperate and semiarid landscapes: Downward approach to water balance 
% analysis. Water Resources Research, 39(2). 
% http://doi.org/10.1029/2001WR000328

theta = [1   , 2000;    % Sb, Maximum soil moisture storage [mm]
         0.05, 0.95;    % phi, Porosity [-]
         0.01, 1.00;    % Sfc, Wilting point as fraction of sb [-]
         0   , 1.00;    % alpha_ss, Subsurface flow constant [1/d]
         1   , 5;       % beta_ss, Subsurface non-linearity constant [-]
         0   , 1.00;    % k_deep, Groundwater recharge constant [d-1]
         0   , 1.00;    % alpha_bf, Groundwater flow constant [d-1]
         1   , 5];      % beta_bf, Groundwater non-linearity constant [-]
      
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'Sb, Maximum soil moisture storage [mm]' ...
                           'phi, Porosity [-]' ...
                           'Sfc, Wilting point as fraction of sb [-]' ...
                           'alpha_ss, Subsurface flow constant [1/d]' ...
                           'beta_ss, Subsurface non-linearity constant [-]' ...
                           'k_deep, Groundwater recharge constant [d-1]' ...
                           'alpha_bf, Groundwater flow constant [d-1]' ...
                           'beta_bf, Groundwater non-linearity constant [-]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)      