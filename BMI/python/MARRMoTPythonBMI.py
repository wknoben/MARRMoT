#! /usr/bin/env python

import os
from oct2py import octave
from BMI import BMI
import numpy as np


class MARRMoTPythonBMI(BMI):
    def __init__(self,MARRMoTLoc=None):
        self._name = "MARRMoTPythonBMI"
        octave_model_root = os.getenv('OCTAVE_MODEL_ROOT', MARRMoTLoc)
        octave.addpath(octave.genpath(octave_model_root))
        octave.eval("model = marrmotBMI_oct()")
               
    def initialize(self, filename):
        commandString = 'model.initialize("' + filename + '")'
        octave.eval(commandString)
        
    def update(self):
        octave.eval("model.update()")
        
    def finalize(self):
        octave.eval("model.finalize()")        

    def get_component_name(self):
        return octave.eval('model.get_component_name()')

    def get_input_var_names(self):
        return octave.eval('model.get_input_var_names()').tolist()[0]

    def get_output_var_names(self):
        return octave.eval('model.get_output_var_names()').tolist()[0]

    def get_var_grid(self, varName):
        commandString = 'model.get_var_grid("' + varName + '")'
        return octave.eval(commandString)

    def get_var_type(self, varName):
        commandString = 'model.get_var_type("' + varName + '")'
        return octave.eval(commandString)

    def get_var_units(self, varName):
        commandString = 'model.get_var_units("' + varName + '")'
        return octave.eval(commandString)

    def get_var_itemsize(self, itemName):
        commandString = 'model.get_var_itemsize("' + itemName + '")'
        return int(octave.eval(commandString))

    def get_var_nbytes(self, varName):
        commandString = 'model.get_var_nbytes("' + varName + '")'
        return int(octave.eval(commandString))

    def get_var_location(self, varName):
        commandString = 'model.get_var_location("' + varName + '")'
        return octave.eval(commandString)

    def get_current_time(self):
        return octave.eval('model.get_current_time()')
 
    def get_start_time(self):
        return octave.eval('model.get_start_time()')
    
    def get_end_time(self):
        return octave.eval('model.get_end_time()')

    def get_time_units(self): 
        return octave.eval('model.get_time_units()')

    def get_time_step(self):
        return octave.eval('model.get_time_step()')

    def get_value(self, varName):
        commandString = 'model.get_value("' + varName + '")'
        value = octave.eval(commandString)
        if type(value) is np.ndarray:
            return value
        else:
            return np.array([value])
            
    def get_value_ptr(self, varName):
        raise NotImplementedError("Reference cannot be transmitted from Octave to Python")
        
    def get_value_at_indices(self, varName, indices):
        # temporary fix since get_value_at_indices doesn't work with indices input and MARRMoT do not have grids
        if str(indices) != '[0, 0]':
            raise Exception('indices are out of bound. The model has no grid, thus indices=[0,0])')    
        else:
            commandString = 'model.get_value("' + varName + '")'
            return octave.eval(commandString)
 
    def set_value(self, varName, src):
        commandString = "model.set_value('" + varName + "'," + np.array2string(src) + ")"
        return octave.eval(commandString)

    def set_value_at_indices(self, varName, indices, src):
        commandString = "model.set_value_at_indices('" + varName + "'," + np.array2string(indices) + "," + np.array2string(src)+ ")"
        return octave.eval(commandString)

    # Grid information
    def get_grid_rank(self, grid_id):
        return octave.eval('model.get_grid_rank(' + str(grid_id) + ')')


    def get_grid_size(self, grid_id):
        return octave.eval('model.get_grid_size(' + str(grid_id) + ')')


    def get_grid_type(self, grid_id):
        return octave.eval('model.get_grid_type(' + str(grid_id) + ')')

    # Uniform rectilinear
    def get_grid_shape(self, grid_id):
        return octave.eval('model.get_grid_shape(' + str(grid_id) + ')').flatten()

    def get_grid_spacing(self, grid_id):
        return octave.eval('model.get_grid_spacing(' + str(grid_id) + ')').flatten()

    def get_grid_origin(self, grid_id):
        return octave.eval('model.get_grid_origin(' + str(grid_id) + ')').flatten()

    # Non-uniform rectilinear, curvilinear
    def get_grid_x(self, grid_id):
        value = octave.eval('model.get_grid_x(' + str(grid_id) + ')')
        # oct2py converts single value vectors to scalars, while bmi expects list
        if type(value) is float:
            return [value]
        return value

    def get_grid_y(self, grid_id):
        value = octave.eval('model.get_grid_y(' + str(grid_id) + ')')
        # oct2py converts single value vectors to scalars, while bmi expects list
        if type(value) is float:
            return [value]
        return value

    def get_grid_z(self, grid_id):
        return octave.eval('model.get_grid_z(' + str(grid_id) + ')')

   #not implemented in MARRMoT (all MARRMoT models have no grid)

    def get_grid_node_count(self, grid_id):
        return "grids not implemented yet"

    def get_grid_edge_count(self, grid_id):
        return "grids not implemented yet"

    def get_grid_face_count(self, grid_id):
        return "grids not implemented yet"

    def get_grid_edge_nodes(self, grid_id):
        return "grids not implemented yet"

    def get_grid_face_nodes(self, grid_id):
        return "grids not implemented yet"

    def get_grid_nodes_per_face(self, grid_id):
        return "grids not implemented yet"
