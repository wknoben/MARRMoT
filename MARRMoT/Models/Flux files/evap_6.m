function [func] = evap_6(~)
%evap_6 Creates function for evaporation: evaporates based on scaled
%current water storage, a wilting point, a constraining factor and limited 
%by potential rate.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
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
%
% WK, 05/10/2018

func = @(p1,p2,S,Smax,Ep,dt) min([p1.*Ep,p1*Ep*S./(p2*Smax),S/dt]);

end

