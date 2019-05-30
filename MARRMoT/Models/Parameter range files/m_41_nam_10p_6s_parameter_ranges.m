function [ theta ] = m_41_nam_10p_6s_parameter_ranges( )
%m_41_nam_10p_6s_parameter_ranges Provides parameter ranges for calibration
%   of the 6-store NAM model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Nielsen, S. A., & Hansen, E. (1973). Numerical simulation of he rainfall-
% runoff process on a daily basis. Nordic Hydrology, (4), 171–190. 
% http://doi.org/https://doi.org/10.2166/nh.1973.0013

theta = [   0, 20;      % cs, Degree-day factor for snowmelt [mm/oC/d]
            0, 1;       % cif, Runoff coefficient for interflow [d-1]
            1, 2000;    % stot, Maximum total soil moisture depth [mm]
            0, 0.99;    % cl1, Lower zone filling threshold for interflow generation [-]
            0.01, 0.99; % f1, Fraction of total soil depth that makes up lstar
            0, 1;       % cof, Runoff coefficient for overland flow [d-1]
            0, 0.99;    % cl2, Lower zone filling threshold for overland flow generation [-]
            0, 1;       % k0, Overland flow routing delay [d-1]
            0, 1;       % k1, Interflow routing delay [d-1]
            0, 1];      % kb, Baseflow routing delay [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'cs, Degree-day factor for snowmelt [mm/oC/d]' ...
                           'cif, Runoff coefficient for interflow [d-1]' ...
                           'stot, Maximum total soil moisture depth [mm]' ...
                           'cl1, Lower zone filling threshold for interflow generation [-]' ...
                           'f1, Fraction of total soil depth that makes up lstar' ...
                           'cof, Runoff coefficient for overland flow [d-1]' ...
                           'cl2, Lower zone filling threshold for overland flow generation [-]' ...
                           'k0, Overland flow routing delay [d-1]' ...
                           'k1, Interflow routing delay [d-1]' ...
                           'kb, Baseflow routing delay [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)        