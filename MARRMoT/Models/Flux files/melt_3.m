function [out] = melt_3(p1,p2,T,S1,S2,St,dt,varargin)
%melt_3 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Glacier melt provided no snow is stored on the ice layer
% Constraints:  f <= S1/dt
% @(Inputs):    p1   - degree-day factor [mm/oC/d]
%               p2   - temperature threshold for snowmelt [oC]
%               T    - current temperature [oC]
%               S1   - current storage in glacier [mm]
%               S2   - current storage in snowpack [mm]
%               St   - storage in S2 threshold below which glacier melt occurs [mm]
%               dt   - time step size [d]
%               varargin(1) - smoothing variable r (default 0.01)
%               varargin(2) - smoothing variable e (default 5.00)

if size(varargin,2) == 0
    out = min(max(p1*(T-p2),0),S1/dt).*smoothThreshold_storage_logistic(S2,St);
elseif size(varargin,2) == 1
    out = min(max(p1*(T-p2),0),S1/dt).*smoothThreshold_storage_logistic(S2,St,varargin(1));
elseif size(varargin,2) == 2
    out = min(max(p1*(T-p2),0),S1/dt).*smoothThreshold_storage_logistic(S2,St,varargin(1),varargin(2));    
end

end

