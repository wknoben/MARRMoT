function [ theta ] = m_07_gr4j_4p_2s_parameter_ranges( )
%m_07_gr4j_cont_4p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the GR4J model (Perrin et al, 2003). 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
%   Perrin, C., Michel, C. and Andreassion, V. (2003). Improvement of a 
%   parsimonious model for streamflow simulation. Journal of Hydrology 279
%   275-289. DOI: 10.1016/S0022-1694(03)00225-7.


theta = [1   , 2000;    % x1 [mm]
        -10  , 15;      % x2 [mm/d]
         1   , 300;     % x3 [mm]
         1   , 15];     % x4 [d]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'x1 [mm]' ...
                           'x2 [mm/d]' ...
                           'x3 [mm]' ...
                           'x4 [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     
     