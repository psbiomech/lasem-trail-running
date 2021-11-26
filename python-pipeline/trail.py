# -*- coding: utf-8 -*-
"""
Main: extract TRAIL test data

@author: Prasanna Sritharan
"""

#import numpy as np
import usersettings as uset
import c3dextract as c3dex
import labsetup as labs
import opensimsetup as osimsetup
import builddatabase as bd



# set the lab
lasem = labs.LabKeyLasemTrail()

# user settings
user = uset.TRAILSettings()

# build output database
traildb = bd.build_database(user)

# extract C3D data and create OpenSim input files
c3dex.c3d_batch_extract(traildb, lasem, 2, 15, user.refmodelfile)

# write OpenSim setup files
#data = osimsetup.write_ground_forces_mot_file(osimkey)
