function [ theta ] = m_46_classic_12p_8s_parameter_ranges( )
%m_46_classic_12p_8s_parameter_ranges Provides parameter ranges for calibration
%   of the 8-store CLASSIC model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Crooks, S. M., & Naden, P. S. (2007). CLASSIC: a semi-distributed 
% rainfall-runoff modelling system. Hydrology and Earth System Sciences, 
% 11(1), 516–531. http://doi.org/10.5194/hess-11-516-2007

theta = [   0, 1;       % fap, Fraction of catchment area that has permeable soils [-]
            0.01, 0.99; % fdp, Fraction of depth of permeable soil that is store Px [-]
            1, 2000;    % dp, Depth of permeable soil [mm]
            0, 1;       % cq, Runoff coefficient for permeable soil [d-1]
            0, 1;       % d1, Fraction of Ps that infiltrates into semi-permeable soil [-]
            0, 1;       % fas, Fraction of (1-fap) that is fas [-]
            0.01, 0.99; % fds, Fraction of depth of semi-permeable soil that is store Sx [-]
            1, 2000;    % ds, Depth of semi-permeable soil [mm]
            0, 1;       % d2, Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]
            0, 1;       % cxq, Quick runoff coefficient for semi-permeable soil [d-1]
            0, 1;       % cxs, Slow runoff coefficient for semi-permeable soil [d-1]
            0, 1];      % cu, Runoff coefficient for impermeable soil [d-1]
        
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'fap, Fraction of catchment area that has permeable soils [-]' ...
                           'fdp, Fraction of depth of permeable soil that is store Px [-]' ...
                           'dp, Depth of permeable soil [mm]' ...
                           'cq, Runoff coefficient for permeable soil [d-1]' ...
                           'd1, Fraction of Ps that infiltrates into semi-permeable soil [-]' ...
                           'fas, Fraction of (1-fap) that is fas [-]' ...
                           'fds, Fraction of depth of semi-permeable soil that is store Sx [-]' ...
                           'ds, Depth of semi-permeable soil [mm]' ...
                           'd2, Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]' ...
                           'cxq, Quick runoff coefficient for semi-permeable soil [d-1]' ...
                           'cxs, Slow runoff coefficient for semi-permeable soil [d-1]' ...
                           'cu, Runoff coefficient for impermeable soil [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)