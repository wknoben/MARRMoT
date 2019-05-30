function [ theta ] = m_11_collie3_6p_2s_parameter_ranges( )
%m_11_collie3_6p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the Collie River Basin 2 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Jothityangkoon, C., M. Sivapalan, and D. Farmer (2001), Process controls
% of water balance variability in a large semi-arid catchment: downward 
% approach to hydrological model development. Journal of Hydrology, 254,
% 174198. doi: 10.1016/S0022-1694(01)00496-6.

theta = [1   , 2000;      % Smax [mm]
         0.05, 0.95;      % fc as fraction of Smax [-] 
         0   , 1 ;        % a, subsurface runoff coefficient [d-1]
         0.05, 0.95;      % M, fraction forest cover [-]
         1   , 5          % b, flow non-linearity [-]
         0,   1];         % lambda, flow distribution [-]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'Smax [mm]' ...
                           'fc as fraction of Smax [-]' ...
                           'a, subsurface runoff coefficient [d-1]' ...
                           'M, fraction forest cover [-]' ...
                           'b, flow non-linearity [-]' ...
                           'lambda, flow distribution [-]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     