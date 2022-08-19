# -*- coding: utf-8 -*-
"""
Build output database for OpenSim projects

@author: Prasanna Sritharan
"""

import glob
import os
import re
import shutil
import pickle as pk



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES



'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
build_database(user, task):
    Build OpenSim output database from user, return metadata dict.
'''
def build_database(user, task):
    
    # parse database, get subject list
    inpath = os.path.join(user.rootpath, user.infolder, user.subjprefix + "*")
    folderlist = glob.glob(inpath)
    subjlist = [os.path.split(f)[1] for f in folderlist]
    
    # build metadata dict
    meta = {}
    outpath = os.path.join(user.rootpath, user.outfolder)
    fnpatobj = re.compile(user.fnpat)   
    for n, subj in enumerate(subjlist):        
        
        # basic info
        meta[subj] = {}
        meta[subj]["subj"] = subj
        meta[subj]["project"] = user.project
        meta[subj]["outpath"] = os.path.join(outpath, subj)
        
        # trial subfolders
        meta[subj]["trials"] = {}
        for group in user.trialgroupfolders:
            
            # parse subfolders
            meta[subj]["trials"][group] = {}
            groupinpath = os.path.join(folderlist[n], group, "*.c3d")
            groupfolderlist = glob.glob(groupinpath)
            triallist = [os.path.splitext(os.path.split(f)[1])[0] for f in groupfolderlist]
            
            # trials (for selected task only)
            for m, trial in enumerate(triallist):
                trialprefix = fnpatobj.fullmatch(trial).group(1)              
                if (trialprefix.casefold() == user.staticprefix.casefold()) or (trialprefix.casefold() in [t.casefold() for t in user.trialprefixes[task.casefold()]]):                                   
                    meta[subj]["trials"][group][trial] = {}                    
                    meta[subj]["trials"][group][trial]["trial"] = trial
                    meta[subj]["trials"][group][trial]["c3dfile"] = trial + ".c3d"
                    meta[subj]["trials"][group][trial]["osim"] = subj + ".osim"
                    meta[subj]["trials"][group][trial]["inpath"] = os.path.split(groupfolderlist[m])[0]
                    meta[subj]["trials"][group][trial]["outpath"] = os.path.join(outpath, subj, group, trial)
                    meta[subj]["trials"][group][trial]["task"] = task
                    meta[subj]["trials"][group][trial]["condition"] = trialprefix.casefold()
                    meta[subj]["trials"][group][trial]["isstatic"] = False
                    meta[subj]["trials"][group][trial]["usedstatic"] = False
                    if trialprefix.casefold() == user.staticprefix.casefold():
                        meta[subj]["trials"][group][trial]["task"] = "static"
                        meta[subj]["trials"][group][trial]["condition"] = "static"
                        meta[subj]["trials"][group][trial]["isstatic"] = True
                        if trial.casefold().endswith(user.staticused.casefold()):
                            meta[subj]["trials"][group][trial]["usedstatic"] = True
                            
    # create subdfolders if required and copy C3D files into output database
    if not os.path.exists(outpath): os.makedirs(outpath)
    for subj in meta:
        if subj == "project": continue
        if not os.path.exists(meta[subj]["outpath"]): os.makedirs(meta[subj]["outpath"])
        for group in meta[subj]["trials"]:
            if not os.path.exists(os.path.join(meta[subj]["outpath"], group)): os.makedirs(os.path.join(meta[subj]["outpath"], group))
            for trial in  meta[subj]["trials"][group]:
                trialoutpath = meta[subj]["trials"][group][trial]["outpath"]
                if not os.path.exists(trialoutpath): os.makedirs(trialoutpath)
                shutil.copy(os.path.join(meta[subj]["trials"][group][trial]["inpath"], meta[subj]["trials"][group][trial]["c3dfile"]), trialoutpath)
                
    # save the metadata dict
    with open(os.path.join(outpath, user.project + ".pkl"),"wb") as fid: pk.dump(meta, fid)
    
    return meta
