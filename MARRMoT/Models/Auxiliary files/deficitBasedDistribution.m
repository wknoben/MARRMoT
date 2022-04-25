function [ f1,f2 ] = deficitBasedDistribution( S1,S1max,S2,S2max )
%DEFICITBASEDDISTRIBUTION Calculates a fractional split for two stores,
%based on the relative deficit in each. Currently used in:
% m_33_sacramento_11p_5s

% Copyright (C) 2018 Wouter J.M. Knoben
% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT).
% MARRMoT is a free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

%%CALCULATE RELATIVE DEFICITS
rd1 = (S1-S1max)/S1max;
rd2 = (S2-S2max)/S2max;

%%CALCULATE FRACTIONAL SPLIT
if rd1+rd2 ~= 0
    
    % Deficit exists and can be used to compute the split
    f1 = rd1/(rd1+rd2);
    f2 = rd2/(rd1+rd2);
else 
    
    % Both deficits are zero, and we revert to distribution based on
    % relative maximum store size
    f1 = S1max/(S1max+S2max);
    f2 = S2max/(S1max+S2max);
end
    
end

