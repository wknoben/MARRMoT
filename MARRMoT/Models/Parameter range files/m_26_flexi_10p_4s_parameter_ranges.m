function [ theta ] = m_26_flexi_10p_4s_parameter_ranges( )
%m_26_flexi_10p_4s_parameter_ranges Provides parameter ranges for calibration
%   of the 4-store Flex-I model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Fenicia, F., McDonnell, J. J., & Savenije, H. H. G. (2008). Learning from
% model improvement: On the contribution of complementary data to process 
% understanding. Water Resources Research, 44(6), 1–13. 
% http://doi.org/10.1029/2007WR006386

theta = [1,2000;        % URmax, Maximum soil moisture storage [mm]
         0, 10;         % beta, Unsaturated zone shape parameter [-]
         0, 1;          % D, Fast/slow runoff distribution parameter [-]
         0, 20;         % PERCmax, Maximum percolation rate [mm/d]
         0.05, 0.95;    % Lp, Wilting point as fraction of s1max [-]
         1, 5;          % Nlagf, Flow delay before fast runoff [d]
         1, 15;         % Nlags, Flow delay before slow runoff [d]
         0, 1;          % Kf, Fast runoff coefficient [d-1]
         0, 1;          % Ks, Slow runoff coefficient [d-1]
         0, 5];         % Imax, Maximum interception storage [mm]
     

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'URmax, Maximum soil moisture storage [mm]' ...
                           'beta, Unsaturated zone shape parameter [-]' ...
                           'D, Fast/slow runoff distribution parameter [-]' ...
                           'PERCmax, Maximum percolation rate [mm/d]' ...
                           'Lp, Wilting point as fraction of s1max [-]' ...
                           'Nlagf, Flow delay before fast runoff [d]' ...
                           'Nlags, Flow delay before slow runoff [d]' ...
                           'Kf, Fast runoff coefficient [d-1]' ...
                           'Ks, Slow runoff coefficient [d-1]' ...
                           'Imax, Maximum interception storage [mm]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     