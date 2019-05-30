function [ theta ] = m_14_topmodel_7p_2s_parameter_ranges( )
%m_14_topmodel_7p_2s_parameter_ranges Provides parameter ranges for calibration
%   of the 2-store TOPMODEL.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Beven, K., Lamb, R., Quinn, P., Romanowicz, R., & Freer, J. (1995). 
% TOPMODEL. In V. P. Singh (Ed.), Computer Models of Watershed Hydrology 
% (pp. 627–668). Baton Rouge: Water Resources Publications, USA.
%
% Beven, K. J., & Kirkby, M. J. (1979). A physically based, variable 
% contributing area model of basin hydrology / Un modèle à base physique 
% de zone d’appel variable de l'hydrologie du bassin versant. Hydrological 
% Sciences Bulletin, 24(1), 43–69. http://doi.org/10.1080/02626667909491834
%
% Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. a., Vrugt, J. a., 
% Gupta, H. V., … Hay, L. E. (2008). Framework for Understanding Structural
% Errors (FUSE): A modular framework to diagnose differences between 
% hydrological models. Water Resources Research, 44(12). 
% http://doi.org/10.1029/2007WR006735
%
% Sivapalan, M., Beven, K., & Wood, E. F. (1987). On hydrologic similarity:
% 2. A scaled model of storm runoff production. Water Resources Research, 
% 23(12), 2266–2278. http://doi.org/10.1029/WR023i012p02266

theta = [1, 2000;    % suzmax, Maximum soil moisture storage in unsatured zone [mm]
         0.05, 0.95; % st, Threshold for flow generation and evap change as fraction of suzmax [-]
         0,1;        % kd, Leakage to saturated zone flow coefficient [mm/d]
         0.1,200;    % q0, Zero deficit base flow speed [mm/d]
         0,1;        % m, Baseflow coefficient [mm-1]
         1, 7.5;     % chi, Gamma distribution parameter [-]
         0.1, 5];    % phi, Gamma distribution parameter [-]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'suzmax, Maximum soil moisture storage in unsatured zone [mm]' ...
                           'st, Threshold for flow generation and evap change as fraction of suzmax [-]' ...
                           'kd, Leakage to saturated zone flow coefficient [mm/d]' ...
                           'q0, Zero deficit base flow speed [mm/d]' ...
                           'm, Baseflow coefficient [mm-1]' ...
                           'chi, Gamma distribution parameter [-]' ...
                           'phi, Gamma distribution parameter [-]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     