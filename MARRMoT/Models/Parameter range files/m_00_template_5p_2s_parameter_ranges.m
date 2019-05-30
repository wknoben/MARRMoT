function [ theta ] = m_00_template_5p_2s_parameter_ranges( )
%m_00_template_5p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store test model, created by W. Knoben in 02-2018.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

theta = [1   , 40;      % Smax [mm]
         0   , 2 ;      % kc, capillary rise [mm/d]
         0   , 3 ;      % kp, percolation rate [mm/d]
         0.5 , 1;       % ks, base flow time parameter [d-1]
         1   , 5];      % time delay of routing scheme [d]
   
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'Smax [mm]' ...
                           'kc, capillary rise [mm/d]' ...
                           'kp, percolation rate [mm/d]' ...
                           'ks, base flow time parameter [d-1]' ...
                           'time delay of routing scheme [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)