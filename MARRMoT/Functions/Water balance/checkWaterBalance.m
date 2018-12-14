function [ out ] = checkWaterBalance( P,fluxes,stores,sini,R )
%checkWaterBalance Calculates a water balance
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Input:
% P         - time series of (precipitation) input              [mm]
% fluxes.Q  - time series of streamflow output                  [mm]
% fluxes.Ea - time series of evaporation                        [mm]
% stores.Sn - structure with time series of n stores            [mm]
% R         - time series of water still in the routing module  [mm]
%
% Note: R is not always present. Input value 0 to indicate this.

% Get variables
numStores = length(sini);
S         = zeros(numStores,1);
Q         = fluxes.Q;
E         = fluxes.Ea;

% Display
disp(['Total P  = ',num2str(sum(P)),' mm.'])
disp(['Total Q  = ',num2str(sum(Q)),' mm.'])
disp(['Total E  = ',num2str(sum(E)),' mm.'])

for s = 1:numStores
    disp(['Delta S',num2str(s),' = ',num2str((stores.(['S',num2str(s)])(end)-sini(s))),' mm.'])
    S(s) = stores.(['S',num2str(s)])(end)-sini(s);
end

out =   sum(P) - ...
        sum(E) - ...
        sum(Q) - ...
        sum(S) - ...
        sum(R);

if length(R) == 1  && R == 0                                               % Routing has placeholder value 0
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a)) - delta S = ',...
          num2str(out),' mm.'])
else                                                                        % Routing is included
    disp(['On route = ',num2str(sum(R)),' mm.'])
    disp(['Water balance = sum(P) - (sum(Q) + sum(E_a)) - delta S - still-being-routed = ',...
          num2str(out),' mm.'])
end
end

