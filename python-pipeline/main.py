# -*- coding: utf-8 -*-
"""
Process and run LASEM TRAIL project data through OpenSim

@author: Prasanna Sritharan
"""

import usersettings as uset
import c3dextract as c3dex
import labsetup as labs
import builddatabase as bd
import pickle as pk
import os


print("\n\n\n")
print("----------------------------------------")
print("TRAIL: DATA PROCESSING & OPENSIM")
print("----------------------------------------")
print("\n")


# %% SET THE LAB

print("Loading lab info...", end="")
lasem = labs.LabKeyLasemTrail()
print("Done.\n")


# %% USER SETTINGS
print("Loading user settings... ", end="")
user = uset.TRAILSettings()
print("Done.\n")


# %% META DATABASE (BUILD NEW OR LOAD EXISTING)

#print("Building new output database... ", end="")
#traildb = bd.build_database("TRAIL", user, "run")
#print("Done.\n")

print("Loading existing output database... ", end="")
dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath,"rb") as fid:
    traildb = pk.load(fid)
print("Done.\n")


# %% EXTRACT C3D AND CREATE OPENSIM DATA FILES

print("Extracting C3D data, creating OpenSim files...\n")
osimkey = c3dex.c3d_batch_process(user, traildb, lasem, 2, 15, -1)
print("\nC3D data extract done.\n")


# %% RUN OPENSIM PIPELINE

# TBD
