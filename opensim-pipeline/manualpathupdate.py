# -*- coding: utf-8 -*-
"""
Manually update paths in output files

@author: Prasanna Sritharan
"""

import pickle as pk
import os
import usersettings as us


# Tasks
task = "hfd"
dataset = "hfd"

# Path roots
oldroot = r"C:\Users\Owner\Documents\data\TRAIL"  # Lenovo
rootpath = r"C:\Users\psrit\Documents\data\TRAIL" # MSI


# Get user settings
print("Loading user settings... ", end="")
user = us.TRAILSettings_HFD()
print("Done.\n")


# ###################################
# META DICT

# Load the meta dict
dbfilepath = os.path.join(rootpath, user.outfolder, task, dataset, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    traildb = pk.load(fid)


# Update paths
for subj in traildb.keys():
    
    # Skip
    if subj=="study": continue

    # Update outpath
    traildb[subj]["outpath"] = os.path.join(rootpath, user.outfolder, task, dataset, subj)
    
    # Update trial inpath and outpath
    for group in traildb[subj]["trials"].keys():
        for trial in traildb[subj]["trials"][group].keys():
            traildb[subj]["trials"][group][trial]["inpath"] = os.path.join(rootpath, user.infolder[0], subj, group)
            traildb[subj]["trials"][group][trial]["outpath"] = os.path.join(rootpath, user.outfolder, task, dataset, subj, group, trial)


# Save updated meta dict
outpath = os.path.join(rootpath, user.outfolder, task, dataset)
with open(os.path.join(outpath, user.project + ".pkl"), "wb") as fid:
    pk.dump(traildb, fid)



# ###################################
# DATA STRUCTURES

# Subject
for subj in traildb.keys():
    
    # skip the study info
    if subj.casefold() == "study": continue 

    # Groups
    for group in traildb[subj]["trials"].keys():
        
        # Trials
        for trial in traildb[subj]["trials"][group].keys():
            
            
            # ###################################
            # TRIALKEY

            # No paths are stored in the TrialKey
            
            
            # ###################################
            # OSIMKEY
            
            # load the OsimKey
            pklpath0 = traildb[subj]["trials"][group][trial]["outpath"]
            pklfile0 = traildb[subj]["trials"][group][trial]["trial"] + "_osimkey.pkl"
            if not os.path.exists(os.path.join(pklpath0, pklfile0)): continue
            with open(os.path.join(pklpath0, pklfile0),"rb") as fid0: 
                osimkey = pk.load(fid0)

            # Update paths
            osimkey.outpath = traildb[subj]["trials"][group][trial]["outpath"]
            
            # Pickle updated OsimKey
            with open(os.path.join(osimkey.outpath, pklfile0),"wb") as h0:
                pk.dump(osimkey, h0)
                
                
            # ###################################
            # OSIMRESULTSKEY
            
            # load the OsimResultsKey
            pklpath1 = traildb[subj]["trials"][group][trial]["outpath"]
            pklfile1 = traildb[subj]["trials"][group][trial]["trial"] + "_opensim_results.pkl"     
            if not os.path.exists(os.path.join(pklpath1, pklfile1)): continue
            with open(os.path.join(pklpath1, pklfile1),"rb") as fid1: 
                osimresultskey = pk.load(fid1)

            # Update paths
            osimresultskey.outpath = traildb[subj]["trials"][group][trial]["outpath"]
            
            # Pickle updated OsimResultsKey
            with open(os.path.join(osimresultskey.outpath, pklfile1),"wb") as h1:
                pk.dump(osimresultskey, h1)          