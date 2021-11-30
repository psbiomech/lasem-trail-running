# -*- coding: utf-8 -*-
"""
Main: extract TRAIL test data

@author: Prasanna Sritharan
"""

#import numpy as np
import usersettings as uset
import c3dextract as c3dex
import labsetup as labs
import builddatabase as bd



# set the lab
lasem = labs.LabKeyLasemTrail()

# user settings
user = uset.TRAILSettings()

# build output database
traildb = bd.build_database(user, "run")

# extract C3D data and create OpenSim input files
osimkey = c3dex.c3d_batch_process(user, traildb, lasem, "run", 2, 15)

# write OpenSim setup files
#data = osimsetup.write_ground_forces_mot_file(osimkey)
