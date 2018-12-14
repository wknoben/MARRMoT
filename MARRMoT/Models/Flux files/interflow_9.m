function [func] = interflow_9(~)
%interflow_9 
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Anonymous function
% ------------------
% Description:  Non-linear interflow if storage exceeds a threshold
% Constraints:  f <= S-p2
%               S-p2 >= 0     prevents numerical issues with complex numbers
% @(Inputs):    p1   - time coefficient [d-1]
%               p2   - storage threshold for flow generation [mm]
%               p3   - exponential scaling parameter [-]    
%               S    - current storage [mm]
%               dt   - time step size [d]
%
% WK, 09/10/2018

func = @(S,p1,p2,p3,dt)  min((S-p2)/dt,(p1.*max(S-p2,0)).^p3);

end

