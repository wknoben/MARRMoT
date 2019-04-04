function [model,class] = createBMIModelInstance()
   % returns the input variable and its class
   global model
   
   model = marrmotBMI_oct()
   class = class(model)