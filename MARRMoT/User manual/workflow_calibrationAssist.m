function [ Qsim ] = workflow_calibrationAssist( model, input_climate, input_s0, parameters, input_solver )
%workflow_calibrationAssist Extracts simulated flow from model output
%structure for a given combination of model, parameter set and inputs.
% 
% Author:   Wouter J.M. Knoben
% Date:     20-12-2018
% Contact:  w.j.m.knoben@bristol.ac.uk
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Run model
out = feval(model,...                % Model name
            input_climate,...        % Climate data
            input_s0,...             % Initial storages
            parameters,...           % Parameter values
            input_solver);           % Solver settings

% Extract Qsim
Qsim = out.Q;
        
end

