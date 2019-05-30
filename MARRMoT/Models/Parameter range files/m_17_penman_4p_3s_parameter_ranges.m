function [ theta ] = m_17_penman_4p_3s_parameter_ranges( )
%m_17_penman_4p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the 3-store Penman model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Penman, H. L. (1950). the Dependence of Transpiration on Weather and Soil
% Conditions. Journal of Soil Science, 1(1), 74–89. 
% http://doi.org/10.1111/j.1365-2389.1950.tb00720.x
%
% Wagener, T., Lees, M. J., & Wheater, H. S. (2002). A toolkit for the 
% development and application of parsimonious hydrological models. In Singh,
% Frevert, & Meyer (Eds.), Mathematical Models of Small Watershed Hydrology
% - Volume 2 (pp. 91–139). Water Resources Publications LLC, USA. 

theta = [
        1, 2000;    % smax, Maximum soil moisture storage [mm]
        0, 1;       % phi, Fraction of direct runoff [-]
        0, 1;       % gam, Evaporation reduction in lower zone [-]
        0, 1];      % k1, Runoff coefficient [d-1]
       
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'smax, Maximum soil moisture storage [mm]' ...
                           'phi, Fraction of direct runoff [-]' ...
                           'gam, Evaporation reduction in lower zone [-]' ...
                           'k1, Runoff coefficient [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)