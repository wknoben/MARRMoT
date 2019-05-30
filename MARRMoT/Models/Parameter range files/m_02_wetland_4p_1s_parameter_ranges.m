function [ theta ] = m_02_wetland_4p_1s_parameter_ranges( )
%m_02_wetland_4p_1s_parameter_ranges Provides parameter ranges for calibration
%   of the Wetland model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Savenije, H. H. G. (2010). “Topography driven conceptual modelling 
% (FLEX-Topo).” Hydrology and Earth System Sciences, 14(12), 2681–2692. 
% https://doi.org/10.5194/hess-14-2681-2010

theta = [0   , 5;      % Dw, interception capacity [mm]
         0   , 10;     % Betaw, soil misture distribution parameter [-]
         1   , 2000;   % Swmax, soil misture depth [mm]
         0   , 1];     % kw, base flow time parameter [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'Dw, interception capacity [mm]' ...
                           'Betaw, soil misture distribution parameter [-]' ...
                           'Swmax, soil misture depth [mm]' ...
                           'kw, base flow time parameter [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     