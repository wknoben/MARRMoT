function [ theta ] = m_42_hycymodel_12p_6s_parameter_ranges( )
%m_42_hycymodel_12p_6s_parameter_ranges Provides parameter ranges for calibration
%   of the 6-store HYCYMODEL.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Fukushima, Y. (1988). A model of river flow forecasting for a small 
% forested mountain catchment. Hydrological Processes, 2(2), 167–185.

theta = [
            0, 1;           % c,    Fraction area that is channel [-]
            0, 5;           % imax, Maximum total interception storage [mm]
            0, 1;           % a,    Fraction stem/trunk interception [-]
            0.01, 0.99;     % fi2,  Fraction of total interception that is trunk/stem interception [mm]
            0, 1;           % kin,  Infiltration runoff coefficient [d-1]
            1, 2000;        % D50,  Soil depth where 50% of area contributes to effective flow [mm]
            0.01, 0.99;     % fd16, Soil depth where 16% of area contributes to effective flow [mm]
            1, 2000;        % sbc, Soil depth where evaporation rate starts to decline [mm]
            0, 1;           % kb, Baseflow runoff coefficient [d-1]
            1, 5;           % pb, Baseflow non-linearity [-]
            0, 1;           % kh, Hillslope runoff coefficient [d-1]
            0, 1];          % kc, Channel runoff coefficient [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'c,    Fraction area that is channel [-]' ...
                           'imax, Maximum total interception storage [mm]' ...
                           'a,    Fraction stem/trunk interception [-]' ...
                           'fi2,  Fraction of total interception that is trunk/stem interception [mm]' ...
                           'kin,  Infiltration runoff coefficient [d-1]' ...
                           'D50,  Soil depth where 50% of area contributes to effective flow [mm]' ...
                           'fd16, Soil depth where 16% of area contributes to effective flow [mm]' ...
                           'sbc, Soil depth where evaporation rate starts to decline [mm]' ...
                           'kb, Baseflow runoff coefficient [d-1]' ...
                           'pb, Baseflow non-linearity [-]' ...
                           'kh, Hillslope runoff coefficient [d-1]' ...
                           'kc, Channel runoff coefficient [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)