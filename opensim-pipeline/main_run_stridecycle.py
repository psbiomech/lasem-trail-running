# -*- coding: utf-8 -*-
"""
Process and run LASEM TRAIL project data through OpenSim

@author: Prasanna Sritharan
"""


import datetime


print("\n\n\n")
print("----------------------------------------")
print("TRAIL: DATA PROCESSING & OPENSIM")
print("----------------------------------------")

# start time stamp
ts0 = datetime.datetime.now();
print("Start: %s" % ts0)

print("----------------------------------------")
print("\n")


# %% SET THE LAB

import labsetup as labs

print("Loading lab info...", end="")
lasem = labs.LabKeyLasemTrail()
print("Done.\n")


# %% USER SETTINGS

import usersettings as uset

print("Loading user settings... ", end="")
user = uset.TRAILSettings_RUN()
print("Done.\n")


# %% GET META DATABASE (BUILD NEW OR LOAD EXISTING)...


# OPTION 1: Build new database

# import builddatabase as bd

# print("Building new output database... ", end="")
# traildb = bd.build_database(user, "run", "run_stridecycle")
# print("Done.\n")


# OPTION 2: Load existing database

import pickle as pk
import os

print("Loading existing output database... ", end="")
dbfilepath = os.path.join(user.rootpath, user.outfolder, "run", "run_stridecycle", user.metadatafile)
with open(dbfilepath, "rb") as fid:
    traildb = pk.load(fid)
print("Done.\n")


# %% EXTRACT C3D AND CREATE OPENSIM DATA FILES

# import c3dextract as c3dex

# print("Extracting C3D data, creating OpenSim files...\n")
# failedfiles = c3dex.c3d_batch_process(user, traildb, lasem, 2, restart="TRAIL483")
# print("\nC3D data extract done.\n")


# %% RUN OPENSIM PIPELINE

# import opensimpipeline as osp

# print("Running OpenSim model scaling: SCALE...\n")
# failedstatic = osp.opensim_pipeline(traildb, user, ["scale"])
# print("\nOpenSim model scaling (SCALE) completed.\n")

# print("Running OpenSim analyses: IK, ID...\n")
# osp.opensim_pipeline(traildb, user, ["ik", "id"])
# print("\nOpenSim analyses (IK, ID) completed.\n")

# print("Running OpenSim analyses: BK...\n")
# osp.opensim_pipeline(traildb, user, ["bk"])
# print("\nOpenSim analyses (BK) completed.\n")

# print("Running OpenSim analyses: SO...\n")
# osp.opensim_pipeline(forcedb, user, ["so"])
# print("\nOpenSim analyses (SO) completed.\n")

# print("Running OpenSim analyses: RRA, CMC...\n")
# osp.opensim_pipeline(forcedb, user, ["rra",  "cmc"])
# print("\nOpenSim analyses (RRA, CMC) completed.\n")

# print("Running OpenSim analyses: JR...\n")
# osp.opensim_pipeline(forcedb, user, ["jr"])
# print("\nOpenSim analyses (JR) completed.\n")


# %% LOAD AND FORMAT RESULTS

import opensimresults as osr

print("Converting OpenSim results to Pickle...\n")
osr.opensim_results_batch_process(traildb, ["ik", "id"], user, 101)
print("\nOpenSim results converted to Pickle.\n")

print("Exporting OpenSim results to CSV...\n")
failedfiles = osr.export_opensim_results(traildb, user, ["ik", "id"], 101)
print("CSV export complete.\n")


# %% ADDITIONAL ANALYSES

# import analyses_workpower as anwp

# print("Running wprk and power analyses...\n")
# failedanalyses = anwp.analyses_batch_process(traildb, user, ["jap", "jaw"], True, 4, 6)
# print("Analyses complete.\n")

# print("Exporting work and power analysis results...\n")
# anwp.export_joint_angular_power(traildb, user, 101)
# anwp.export_joint_angular_work(traildb, user)
# print("Analyses results export complete.\n")


# %% END

print("\n")
print("----------------------------------------")

# end time stamp
ts1 = datetime.datetime.now();
print("End: %s" % ts1)
datetime.datetime.now()

print("----------------------------------------")
print("\n")


