function [out] = evap_6(p1,p2,S,Smax,Ep,dt)
%evap_6 evaporation based on scaled current water storage, a wilting point,
%a constraining factor and limited by potential rate.

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Transpiration from vegetation at the potential rate if 
%               storage is above a wilting point and scaled by relative 
%               storage if not
% Constraints:  Ea <= Ep
%               Ea <= S/dt
% @(Inputs):    p1   - fraction vegetated area [-]
%               p2   - wilting point as fraction of Smax
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               Ep   - potential evapotranspiration rate [mm/d]
%               dt   - time step size [d]

out = min([p1.*Ep,p1*Ep*S./(p2*Smax),S/dt]);

end

