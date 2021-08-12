function [ uh ] = update_uh(uh, flux_in)
% UPDATE_UH calculates new still-to-flow values of a unit hydrograph at the
% based on a flux routed through it at this timestep.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% In:
% uh        - unit hydrograph    [nx2]
%               uh's first row contains coeficients to splut flow at each
%               of n timesteps forward, the second row contains
%               still-to-flow values.
% flux_in   - input flux         [1x1]
%
% Out:
% uh        - update unit hydrograph [nx2]
%               UPDATE_UH does not change the first row, i.e. the
%               coefficients, but only the still-to-flow values.
%

uh(2,:)   = (uh(1,:) .* flux_in) + uh(2,:);
uh(2,:)   = circshift(uh(2,:),-1);
uh(2,end) = 0;

end