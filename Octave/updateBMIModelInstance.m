function [model,class] = updateBMIModelInstance()
   % returns the input variable and its class
   global model
   
   model.update()
   class = class(model)