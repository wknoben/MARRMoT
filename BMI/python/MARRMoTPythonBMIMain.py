#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 22:00:23 2019

@author: rwhut
"""

import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import os

os.environ['OCTAVE_EXECUTABLE'] = '/Applications/Octave-4.4.1.app/Contents/Resources/usr/bin/octave'
import MARRMoTPythonBMI


#make instance of BMILorenz
model = MARRMoTPythonBMI.MARRMoTPythonBMI('/Users/rwhut/Documents/eWaterCycle/repos/MARRMoT')

model.initialize('/Users/rwhut/Documents/eWaterCycle/repos/MARRMoT/BMI/Config/BMI_testcase_m01_BuffaloRiver_TN_USA')

discharge = []


while model.get_current_time() < 356:
    model.update()
    flux = model.get_value('flux_out')
    discharge.append(flux['Q'])
    
plt.plot(discharge)
