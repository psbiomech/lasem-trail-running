# -*- coding: utf-8 -*-
"""
Main: extract TRAIL test data

@author: Prasanna Sritharan
"""

#import numpy as np
import c3dextract as c3dex
import labsetup as labs
import opensimsetup as osimsetup


# set the lab
lasem = labs.LabKeyLasemTrail()

# return C3D data
c3dpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data\Events\TRAIL_071\Baseline"
c3dfile = "TRAIL_071_EP_01.c3d"
ref_model = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\reference_model\Rajagopal2015.osim"
c3dkey, trialkey, osimkey = c3dex.c3d_extract(c3dfile, c3dpath, lasem, 2, 15, ref_model)

# write OpenSim setup files
data = osimsetup.write_ground_forces_mot_file(osimkey)
