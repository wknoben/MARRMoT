function [res] = checkTimeseriesSize(s1,s2)
%checkTimeseriesSize Checks if two time series are of equal size.%
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
%   Compares the size of two time series, expects vectors as input. 
%   Possible outputs:
%   1 - time series are equal
%   2 - equal size, but different orientations
%   0 - not equal size

res = 0;

if size(s1) == size(s2)
    res = 1;
elseif size(s1,1) == size(s2,2) && size(s1,2) == size(s2,1)
    res = 2;
end


end

