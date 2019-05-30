function [ theta ] = m_25_tcm_6p_4s_parameter_ranges( )
%m_25_tcm_6p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store Thames Catchment model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Moore, R. J., & Bell, V. A. (2001). Comparison of rainfall-runoff models 
% for flood forecasting. Part 1: Literature review of models. Bristol: 
% Environment Agency.

theta = [   0, 1;     % phi, Fraction preferential recharge [-]
            1, 2000;  % rc, Maximum soil moisture depth [mm]
            0, 1;     % gam, Fraction of Ep reduction with depth [-]
            0, 1;     % k1, Runoff coefficient [d-1]
            0, 1;     % fa, Fraction of mean(P) that forms abstraction rate [mm/d]
            0, 1];    % k2, Runoff coefficient [mm-1 d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'phi, Fraction preferential recharge [-]' ...
                           'rc, Maximum soil moisture depth [mm]' ...
                           'gam, Fraction of Ep reduction with depth [-]' ...
                           'k1, Runoff coefficient [d-1]' ...
                           'fa, Fraction of mean(P) that forms abstraction rate [mm/d]' ...
                           'k2, Runoff coefficient [mm-1 d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)