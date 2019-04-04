#!python
"""
wflow_prepare
=============
A script to run Marrmot octave using oct2py. Still in a very initial development phase.
Usage::
    
.. todo:
    # Currently oct2py is not able to read the methods from the class. May not be possible...
"""

import numpy as np
import scipy
from oct2py import octave

"""
set path o config file
"""
path = '/home/yifat/MARRMoT/BMI/Config/BMI_testcase_m01_BuffaloRiver_TN_USA' 

"""
Call the BMI implementation in octave: marrmotBMI_oct.m
"""
obj_marrmot_m01 = octave.marrmotBMI_oct()

"""
Here we recieve an error that the object has no attribute initialize:
"""
obj_marrmot_m01.initialize(path)

"""
dir(obj_marrmot_m01) will return only the properties, does not show the implemented methods
"""
