function [x,class] = update()
   % returns the input variable and its class
   global state 
   state = state + 1
   x = state
   class = class(x)