classdef BMI < hgsetget
    % a superclass for models that want to implement the Basic Model
    % Interface (BMI, see http://csdms.colorado.edu/wiki/BMI_Description).
    % to construct a model, use: classdef myModel < BMI.
    % in the constructor of your class, the properties defined below should
    % be given values. Note that the name, var_units, input_var_units and
    % output_var_units need all be map containers. See the AR1BMI.m file
    % for an example.
    
    
    %These properties are shared by all BMI models. they need to be set in
    %the constructor of the model class.
    properties
    end
    
    
    %the following is the heart of the BMI class. Every BMI model needs to
    %implement these methods. For some methods a standard implementation is
    %provided (example: update_until). These can be overwritten by
    %implementing them in the model class file. See the AR1class.m for an
    %example. Methods that are not implemented will throw an error msg
    %declaring they are not implemented.
    %See the BMI reference for more details.
    methods
        function initialize(obj,fileName)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function update(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        
        
        function update_until(obj,time)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_var_units(obj,long_var_name)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        
        function output = get_var_type(obj,long_var_name)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_var_rank(obj,long_var_name)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_value(obj,long_var_name)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_value_at_indices(obj,long_var_name, inds)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        
        function set_value(obj,long_var_name, src)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        
        function set_value_at_indices(obj,long_var_name, inds, src)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        
        function output = get_component_name(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_input_var_names(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_output_var_names(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_start_time(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_end_time(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
        function output = get_current_time(obj)
            error(['This method is not defined in class ' mfilename('class')]);
        end
    end
end

