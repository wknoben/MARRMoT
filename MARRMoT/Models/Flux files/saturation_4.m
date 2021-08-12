function [out] = saturation_4(S,Smax,In)
%saturation_4 

% Copyright (C) 2019, 2021 Wouter J.M. Knoben, Luca Trotter
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% Flux function
% ------------------
% Description:  Saturation excess from a store with different degrees of 
%               saturation (quadratic variant)
% Constraints:  f >= 0
% @(Inputs):    S    - current storage [mm]
%               Smax - maximum storage [mm]
%               In   - incoming flux [mm/d]

out = max(0,(1-(S/Smax).^2).*In);

end

