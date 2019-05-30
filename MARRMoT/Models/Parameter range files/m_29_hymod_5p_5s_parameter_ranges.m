function [ theta ] = m_29_hymod_5p_5s_parameter_ranges( )
%m_29_hymod_5p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the 5-store HyMOD model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Boyle, D. P. (2001). Multicriteria calibration of hydrologic models. 
% University of Arizona. Retrieved from http://hdl.handle.net/10150/290657
%
% Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, Hoshin, 
% V., & Sorooshian, S. (2001). A framework for development and application 
% of hydrological models. Hydrology and Earth System Sciences, 5, 13–26.

theta = [1  , 2000;     % Smax, Maximum soil moisture storage [mm], 
         0  , 10;       % b, Soil depth distribution parameter [-]
         0  , 1;        % a, Runoff distribution fraction [-]
         0  , 1;        % kf, base flow time parameter [d-1]
         0  , 1];       % ks, base flow time parameter [d-1]
     
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'Smax, Maximum soil moisture storage [mm], ' ...
                           'b, Soil depth distribution parameter [-]' ...
                           'a, Runoff distribution fraction [-]' ...
                           'kf, base flow time parameter [d-1]' ...
                           'ks, base flow time parameter [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)