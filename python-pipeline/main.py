# -*- coding: utf-8 -*-
"""
Process and run LASEM TRAIL project data through OpenSim

@author: Prasanna Sritharan
"""

import c3dextract as c3dex
import builddatabase as bd
import usersettings as uset
import labsetup as labs
import opensimpipeline as osp
import opensimresults as osr
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

print("Building new output database... ", end="")
traildb = bd.build_database("TRAIL", user, "run_stridecycle")
print("Done.\n")

# print("Loading existing output database... ", end="")
# dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
# with open(dbfilepath,"rb") as fid:
#     traildb = pk.load(fid)
# print("Done.\n")


# %% EXTRACT C3D AND CREATE OPENSIM DATA FILES

print("Extracting C3D data, creating OpenSim files...\n")
c3dex.c3d_batch_process(user, traildb, lasem, 2, -1)
print("\nC3D data extract done.\n")


# %% RUN OPENSIM PIPELINE

print("Running OpenSim basic analyses (IK, ID)...\n")
osp.opensim_pipeline(traildb, user, ["scale", "ik", "id"])
print("\nOpenSim pipeline completed.\n")

# print("Running additional OpenSim analyses (SO)...\n")
# osp.opensim_pipeline(traildb, user, ["so"])
# print("\nOpenSim pipeline completed.\n")

# print("Running additional OpenSim analyses (RRA, CMC)...\n")
# osp.opensim_pipeline(traildb, user, ["rra", "cmc"])
# print("\nOpenSim pipeline completed.\n")

print("Running additional OpenSim analyses (SO + CMC)...\n")
osp.opensim_pipeline(traildb, user, ["so", "rra", "cmc"])
print("\nOpenSim pipeline completed.\n")

# print("Running custom OpenSim analyses series...\n")
# osp.opensim_pipeline(traildb, user, ["scale", "rra", "cmc"])
# print("\nOpenSim pipeline completed.\n")


# %% LOAD AND FORMAT RESULTS

print("Converting OpenSim results to Pickle...\n")
osr.opensim_results_batch_process(traildb, ["scale", "ik", "id", "so"], 101)
print("\nOpenSim results converted to Pickle.\n")


# %% COLLATE RESULTS FOR RSTATS ANALYSIS

print("Exporting OpenSim results to CSV...\n")
osr.export_opensim_results(traildb, user, ["scale", "ik", "id", "so"])
print("CSV export complete.\n")



