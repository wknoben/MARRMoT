function [func] = evap_1(~)
%evap_1 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Evaporation at the potential rate
% Constraints:  f <= S/dt
% @(Inputs):    S    - current storage [mm]
%               Ep   - potential evaporation rate [mm/d]
%               dt   - time step size
%
% WK, 05/10/2018

func = @(S,Ep,dt) min(S/dt,Ep);

end

