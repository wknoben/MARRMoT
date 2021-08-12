function [ flux_out ] = route(flux_in, uh)
% ROUTE calculates the output of a unit hydrograph at the current timestep
% after routing a flux through it.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% flux_in   - input flux         [1x1]
% uh        - unit hydrograph    [nx2]
%               uh's first row contains coeficients to splut flow at each
%               of n timesteps forward, the second row contains
%               still-to-flow values.
%
% Out:
% flux_out  - flux routed through the uh at this step [1x1]
%

flux_out = uh(1,1) * flux_in + uh(2,1);

end