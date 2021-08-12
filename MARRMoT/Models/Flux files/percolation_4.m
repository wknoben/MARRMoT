function [out] = percolation_4(p1,p2,p3,p4,p5,S,Smax,dt)
%percolation_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Demand-based percolation scaled by available moisture
% Constraints:  f <= S/dt
%               f >= 0          prevents erratic numerical behaviour
% @(Inputs):    p1   - base percolation rate [mm/d]
%               p2   - percolation rate increase due moisture deficiencies [mm/d]
%               p3   - non-linearity parameter [-]
%               p4   - summed deficiency across all model stores [mm]
%               p5   - summed capacity of model stores [mm]
%               S    - current storage in the supplying store [mm]
%               Smax - maximum storage in the supplying store [mm]
%               dt   - time step size [d]

% Note: for certain extreme parameter values (very small stores, highly
% non-linear p3) and small computational errors that lead to small negative
% S values, this function behaves erratically. The max(0,S/Smax) part
% prevents this. Similarly, the first max(0,...) part prevents negative
% percolation demands as a result of small numerical errors.

out = max(0,min(S/dt,max(0,S./Smax).*(p1.*(1+p2.*(p4./p5).^(1+p3)))));

end

