function [x, class] = init(y)
   % returns the input variable and its class
   global state = y
   x = state
   class = class(x)