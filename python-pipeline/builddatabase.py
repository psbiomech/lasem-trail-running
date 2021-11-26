# -*- coding: utf-8 -*-
"""
Build output database for OpenSim projects

@author: Prasanna Sritharan
"""

import glob
import os
import shutil
import opensimsetup as osimsetup


'''
build_database(usersettings):
    Build OpenSim output database from UserSettings, return metadata dict.
'''
def build_database(usersettings):
    
    # parse database, get subject list
    inpath = os.path.join(usersettings.rootpath, usersettings.infolder, usersettings.subjprefix + "*")
    folderlist = glob.glob(inpath)
    subjlist = [os.path.split(f)[1] for f in folderlist]
    
    # build metadata dict
    meta = {}
    outpath = os.path.join(usersettings.rootpath, usersettings.outfolder)
    for n, subj in enumerate(subjlist):
        
        # basic info
        meta[subj] = {}
        meta[subj]["subj"] = subj
        meta[subj]["outpath"] = os.path.join(outpath, subj)
        
        # trial subfolders
        meta[subj]["trials"] = {}
        for group in usersettings.trialgroupfolders:
            
            # parse subfolders
            meta[subj]["trials"][group] = {}
            groupinpath = os.path.join(folderlist[n], group, "*.c3d")
            groupfolderlist = glob.glob(groupinpath)
            triallist = [os.path.splitext(os.path.split(f)[1])[0] for f in groupfolderlist]
            
            # trials
            for m, trial in enumerate(triallist):
                meta[subj]["trials"][group][trial] = {}
                meta[subj]["trials"][group][trial]["trial"] = trial
                meta[subj]["trials"][group][trial]["c3dfile"] = trial + ".c3d"
                meta[subj]["trials"][group][trial]["osim"] = subj + ".osim"
                meta[subj]["trials"][group][trial]["inpath"] = os.path.split(groupfolderlist[m])[0]
                meta[subj]["trials"][group][trial]["outpath"] = os.path.join(outpath, subj, group, trial)

    
    # create subdfolders if required and copy C3D files into output database
    if not os.path.exists(outpath): os.makedirs(outpath)
    for subj in meta:
        if not os.path.exists(meta[subj]["outpath"]): os.makedirs(meta[subj]["outpath"])
        for group in meta[subj]["trials"]:
            if not os.path.exists(os.path.join(meta[subj]["outpath"], group)): os.makedirs(os.path.join(meta[subj]["outpath"], group))
            for trial in  meta[subj]["trials"][group]:
                trialoutpath = meta[subj]["trials"][group][trial]["outpath"]
                if not os.path.exists(trialoutpath): os.makedirs(trialoutpath)
                shutil.copy(os.path.join(meta[subj]["trials"][group][trial]["inpath"], meta[subj]["trials"][group][trial]["c3dfile"]), trialoutpath)
                
    return meta
