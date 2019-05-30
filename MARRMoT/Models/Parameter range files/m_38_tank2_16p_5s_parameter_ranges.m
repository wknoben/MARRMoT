function [ theta ] = m_38_tank2_16p_5s_parameter_ranges( )
%m_38_tank2_16p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the 5-store Tank Model - SMA.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Sugawara, M. (1995). Tank model. In V. P. Singh (Ed.), Computer models of
% watershed hydrology (pp. 165–214). Water Resources Publications, USA.

theta = [   0, 1;           % a0, Time parameter for drainage 1>2 [d-1]
            0, 1;           % b0, Time parameter for drainage 2>3 [d-1]
            0, 1;           % c0, Time parameter for drainage 3>4 [d-1]
            0, 1;           % a1, Time parameter for surface runoff 1 [d-1]
            0, 1;           % fa, Fraction of a1 that is a2 [-]
            0, 1;           % fb, Fraction of a2 that is b1 [-]   
            0, 1;           % fc, Fraction of b1 that is c1 [-]
            0, 1;           % fd, Fraction of c1 that is d1 [-]
            1, 2000;        % st, Maximum soil depth (sum of runoff thresholds) [mm]
            0.01, 0.99;     % f2, Fraction of st that consitutes threshold t2 [-]
            0.01, 0.99;     % f1, Fraction of st-t2 that is added to t2 to find threshold 1 [-] (ensures t1 > t2)
            0.01, 0.99;     % f3, Fraction of st-t1-t2 that consitutes threshold 3 [-]
            0, 4;           % k1, Base rate of capillary rise [mm/d]
            0, 4;           % k2, Base rate of soil moisture exchange [mm/d]
            0.01, 0.99;     % z1, Fraction Stot that is sm1 [-]
            0.01, 0.99;];   % z2, Fraction of Stot-sm1 that is sm2 [-]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'a0, Time parameter for drainage 1>2 [d-1]' ...
                           'b0, Time parameter for drainage 2>3 [d-1]' ...
                           'c0, Time parameter for drainage 3>4 [d-1]' ...
                           'a1, Time parameter for surface runoff 1 [d-1]' ...
                           'fa, Fraction of a1 that is a2 [-]' ...
                           'fb, Fraction of a2 that is b1 [-]   ' ...
                           'fc, Fraction of b1 that is c1 [-]' ...
                           'fd, Fraction of c1 that is d1 [-]' ...
                           'st, Maximum soil depth (sum of runoff thresholds) [mm]' ...
                           'f2, Fraction of st that consitutes threshold t2 [-]' ...
                           'f1, Fraction of st-t2 that is added to t2 to find threshold 1 [-] (ensures t1 > t2)' ...
                           'f3, Fraction of st-t1-t2 that consitutes threshold 3 [-]' ...
                           'k1, Base rate of capillary rise [mm/d]' ...
                           'k2, Base rate of soil moisture exchange [mm/d]' ...
                           'z1, Fraction Stot that is sm1 [-]' ...
                           'z2, Fraction of Stot-sm1 that is sm2 [-]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)        