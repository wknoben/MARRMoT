function [ theta ] = m_28_xinanjiang_12p_4s_parameter_ranges( )
%m_28_xinanjiang_12p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store Xinanjiang model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Jayawardena, A. W., & Zhou, M. C. (2000). A modified spatial soil moisture
% storage capacity distribution curve for the Xinanjiang model. Journal of 
% Hydrology, 227(1-4), 93–113. http://doi.org/10.1016/S0022-1694(99)00173-0
% 
% Zhao, R.-J. (1992). The Xinanjiang model applied in China. Journal of 
% Hydrology, 135(1-4), 371–381. http://doi.org/10.1016/0022-1694(92)90096-E

theta = [   0, 1;           % aim,  Fraction impervious area [-]
            -0.49, 0.49;    % a,    Tension water distribution inflection parameter [-]
            0, 10;          % b,    Tension water distribution shape parameter [-]
            1, 2000;        % stot, Total soil moisture storage (W+S) [mm]
            0.01,0.99;      % fwm,  Fraction of Stot that is Wmax [-]
            0.01,0.99;      % flm,  Fraction of wmax that is LM [-]
            0.01, 0.99;     % c,    Fraction of LM for second evaporation change [-]
            0, 10;          % ex,   Free water distribution shape parameter [-]
            0, 1;           % ki,   Free water interflow parameter [d-1]
            0, 1;           % kg,   Free water groundwater parameter [d-1]
            0, 1;           % ci,   Interflow time coefficient [d-1]
            0, 1];          % cg,   Baseflow time coefficient [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'aim,  Fraction impervious area [-]' ...
                           'a,    Tension water distribution inflection parameter [-]' ...
                           'b,    Tension water distribution shape parameter [-]' ...
                           'stot, Total soil moisture storage (W+S) [mm]' ...
                           'fwm,  Fraction of Stot that is Wmax [-]' ...
                           'flm,  Fraction of wmax that is LM [-]' ...
                           'c,    Fraction of LM for second evaporation change [-]' ...
                           'ex,   Free water distribution shape parameter [-]' ...
                           'ki,   Free water interflow parameter [d-1]' ...
                           'kg,   Free water groundwater parameter [d-1]' ...
                           'ci,   Interflow time coefficient [d-1]' ...
                           'cg,   Baseflow time coefficient [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)
