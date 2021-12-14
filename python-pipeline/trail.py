# -*- coding: utf-8 -*-
"""
Process and run LASEM TRAIL project data through OpenSim

@author: Prasanna Sritharan
"""

import usersettings as uset
import c3dextract as c3dex
import labsetup as labs
import builddatabase as bd


print("\n\n\n")
print("----------------------------------------")
print("TRAIL: DATA PROCESSING & OPENSIM")
print("----------------------------------------")
print("\n")

# set the lab
print("Loading lab info...", end="")
lasem = labs.LabKeyLasemTrail()
print("Done.\n")

# user settings
print("Loading user settings... ", end="")
user = uset.TRAILSettings()
print("Done.\n")

# build output database
print("Building output database... ", end="")
traildb = bd.build_database("TRAIL", user, "run")
print("Done.\n")

# extract C3D data and create OpenSim input files
print("Extracting C3D data, creating OpenSim files...\n")
osimkey = c3dex.c3d_batch_process(user, traildb, lasem, 2, 15, -1)
print("\nC3D data extract done.\n")
