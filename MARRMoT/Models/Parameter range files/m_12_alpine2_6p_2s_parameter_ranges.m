function [ theta ] = m_12_alpine2_6p_2s_parameter_ranges( )
%m_12_alpine2_6p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store Alpine model v1.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Eder, G., Sivapalan, M., & Nachtnebel, H. P. (2003). Modelling water 
% balances in an Alpine catchment through exploitation of emergent 
% properties over changing time scales. Hydrological Processes, 17(11), 
% 2125–2149. http://doi.org/10.1002/hyp.1325

theta = [-3  , 5        % tt [celsius]
         0   , 20;      % degree-day-factor [mm/degree celsius/d]
         1   , 2000;    % Smax [mm]
         0.05, 0.95;    % Field capacity as fraction of Smax [-]
         0   , 1;       % time delay of interflow [d-1]
         0   , 1];      % time delay of baseflow [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'tt [celsius]' ...
                           'degree-day-factor [mm/degree celsius/d]' ...
                           'Smax [mm]' ...
                           'Field capacity as fraction of Smax [-]' ...
                           'time delay of interflow [d-1]' ...
                           'time delay of baseflow [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     