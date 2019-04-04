#! /usr/bin/env python

from oct2py import octave
from BMI import BMI
import numpy as np




class MARRMoTPythonBMI(BMI):
    def __init__(self,MARRMoTLoc):
        self._name = "MARRMoTPythonBMI"
        octave.addpath(octave.genpath(MARRMoTLoc))
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
        return octave.eval('model.get_input_var_names()')

    def get_output_var_names(self):
        return octave.eval('model.get_output_var_names()')

    def get_var_grid(self, gridType):
        commandString = 'model.get_var_grid("' + gridType + '")'
        return octave.eval(commandString)

    def get_var_type(self, varName):
        commandString = 'model.get_var_type("' + varName + '")'
        return octave.eval(commandString)

    def get_var_units(self, varName):
        commandString = 'model.get_var_units("' + varName + '")'
        return octave.eval(commandString)

    def get_var_itemsize(self, itemName):
        commandString = 'model.get_var_itemsize("' + itemName + '")'
        return octave.eval(commandString)

    def get_var_nbytes(self, varName):
        commandString = 'model.get_var_nbytes("' + varName + '")'
        return octave.eval(commandString)

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
        return octave.eval(commandString)
    
    def get_value_ptr(self, varName):
        commandString = 'model.get_value_ptr("' + varName + '")'
        return octave.eval(commandString)

    def get_value_at_indices(self):
        return "not implemented yet"

    def set_value(self, varName, src):
        commandString = "model.set_value(" + varName + "," + np.array2string(src) + ")"
        return octave.eval(commandString)

    def set_value_at_indices(self, varName, indices, src):
        commandString = "model.set_value_at_indices(" + varName + "," + np.array2string(indices) + "," + np.array2string(src)+ ")"
        return octave.eval(commandString)

    # Grid information
    def get_grid_rank(self):
        return "grids not implemented yet"

    def get_grid_size(self):
        return "grids not implemented yet"

    def get_grid_type(self):
        return "grids not implemented yet"

    # Uniform rectilinear
    def get_grid_shape(self):
        return "grids not implemented yet"

    def get_grid_spacing(self):
        return "grids not implemented yet"

    def get_grid_origin(self):
        return "grids not implemented yet"

    # Non-uniform rectilinear, curvilinear
    def get_grid_x(self):
        return "grids not implemented yet"

    def get_grid_y(self):
        return "grids not implemented yet"

    def get_grid_z(self):
        return "grids not implemented yet"

    def get_grid_node_count(self):
        return "grids not implemented yet"

    def get_grid_edge_count(self):
        return "grids not implemented yet"

    def get_grid_face_count(self):
        return "grids not implemented yet"

    def get_grid_edge_nodes(self):
        return "grids not implemented yet"

    def get_grid_face_nodes(self):
        return "grids not implemented yet"

    def get_grid_nodes_per_face(self):
        return "grids not implemented yet"
