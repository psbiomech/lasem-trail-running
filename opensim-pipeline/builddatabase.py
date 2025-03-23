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
import pandas as pd



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
build_database(user, task, dataset, emgosubcohort):
    Build OpenSim output database from user, return metadata dict.
'''
def build_database(user, task, dataset, emgsubcohort=False):
 
    # meta dict
    meta = {}
    meta["study"] = {}
    meta["study"]["task"] = task
    meta["study"]["dataset"] = dataset
 
    # If EMG subcohort only, then load list of subjects, also get group, sex and 
    # symptomatic limb
    if emgsubcohort:
        print("Using EMG subcohort only...")
        emginfo = pd.read_excel(os.path.join(user.rootpath, user.emglistfile), usecols="B, D, E, G")
        emglist = ["TRAIL" + s[-3:] for s in emginfo["ID"].dropna()]
        sexlist = emginfo["Sex"].dropna().to_list()
        typelist = emginfo["Group"].dropna().to_list()
        kneelist = emginfo["knee_reference"].dropna().to_list()
        
 
    
    # input folders
    for infolder in user.infolder:
 
        # parse database, get subject list
        inpath = os.path.join(user.rootpath, infolder, user.subjprefix + "*")
        folderlist0 = glob.glob(inpath, recursive = True)
        subjlist0 = [os.path.split(f)[1] for f in folderlist0]
        
        # If EMG subcohort then cut down to only required subjects
        if emgsubcohort:
            subjidx = [xn for xn, x in enumerate(subjlist0) if x in emglist]
            subjlist = [subjlist0[x] for x in subjidx]    
            folderlist = [folderlist0[x] for x in subjidx]  
        else:
            subjlist = subjlist0
            folderlist = folderlist0
        
        # build metadata dict
        outpath = os.path.join(user.rootpath, user.outfolder, task, dataset)
        fnpatobj = re.compile(user.fnpat) 
        for n, subj in enumerate(subjlist):        
            
            # basic info
            meta[subj] = {}
            meta[subj]["subj"] = subj
            meta[subj]["project"] = user.project
            meta[subj]["outpath"] = os.path.join(outpath, subj)
            meta[subj]["type"] = typelist[n].lower()
            meta[subj]["sex"] = sexlist[n].lower()
            meta[subj]["knee"] = kneelist[n].lower()
            
            # trial subfolders
            meta[subj]["trials"] = {}
            for group in user.trialgroupfolders:
                
                # parse subfolders
                meta[subj]["trials"][group] = {}
                groupinpath = os.path.join(folderlist[n], group, "*.c3d")
                groupfolderlist = glob.glob(groupinpath)
                triallist = [os.path.splitext(os.path.split(f)[1])[0] for f in groupfolderlist]
                
                # skip groups that don't have a least one static trial and at
                # least one required dynamic file
                hasstatic = any([user.staticprefix.casefold() in t.casefold() for t in triallist])
                hasdynamic = any([any([c.casefold() in t.casefold() for c in user.trialprefixes[task][dataset]]) for t in triallist])
                if not(hasstatic) or not(hasdynamic): continue
                
                # trials (for selected dataset only)
                for m, trial in enumerate(triallist):            
                    
                    trial = trial.upper()
                    
                    # file name tokens, skip those that don't match regex
                    trialtoks = fnpatobj.fullmatch(trial)
                    if trialtoks is None: continue
                    
                    # build meta data dict
                    trialprefix = trialtoks.group(user.tasktoknum)                    
                    if (trialprefix.casefold() == user.staticprefix.casefold()) or (trialprefix.casefold() in [t.casefold() for t in user.trialprefixes[task.casefold()][dataset.casefold()]]):                                   
                        meta[subj]["trials"][group][trial] = {}                    
                        meta[subj]["trials"][group][trial]["trial"] = trial
                        meta[subj]["trials"][group][trial]["c3dfile"] = trial + ".c3d"
                        meta[subj]["trials"][group][trial]["osim"] = subj.upper() + ".osim"
                        meta[subj]["trials"][group][trial]["inpath"] = os.path.split(groupfolderlist[m])[0]
                        meta[subj]["trials"][group][trial]["outpath"] = os.path.join(outpath, subj, group, trial)
                        meta[subj]["trials"][group][trial]["task"] = task                        
                        meta[subj]["trials"][group][trial]["dataset"] = dataset
                        meta[subj]["trials"][group][trial]["condition"] = trialprefix.casefold()
                        meta[subj]["trials"][group][trial]["isstatic"] = False
                        meta[subj]["trials"][group][trial]["usedstatic"] = False
                        if trialprefix.casefold() == user.staticprefix.casefold():
                            meta[subj]["trials"][group][trial]["dataset"] = "static"
                            meta[subj]["trials"][group][trial]["condition"] = "static"
                            meta[subj]["trials"][group][trial]["isstatic"] = True
                                                   
                # determine which file to use as static trial in OpenSim, in
                # most cases use the static file set in the user settings, if
                # not found, then use the first available                        
                hasusedstatic = any([meta[subj]["trials"][group][t]["usedstatic"] for t in meta[subj]["trials"][group].keys()])
                if not hasusedstatic:
                    for trial in meta[subj]["trials"][group]:             
                        if trial.casefold().endswith(user.staticused.casefold()):
                            meta[subj]["trials"][group][trial]["usedstatic"] = True                            
                            break
                        elif user.staticprefix.casefold() in trial.casefold():
                            meta[subj]["trials"][group][trial]["usedstatic"] = True
                            break
                    
        # clean up, remove empty subject meta dict keys or those without at
        # least one static and one dynamic trial
        for subj in subjlist:
            if subj.casefold() == "study": continue
            for group in user.trialgroupfolders:
                hasstatic = any([meta[subj]["trials"][group][t]["isstatic"] for t in meta[subj]["trials"][group].keys()])
                hasdynamic = any([not meta[subj]["trials"][group][t]["isstatic"] for t in meta[subj]["trials"][group].keys()])
                if not (hasstatic and hasdynamic):            
                    meta[subj]["trials"].pop(group)
            if not meta[subj]["trials"]:
                meta.pop(subj)
                
        # create subdfolders if required, copy C3D files into output database
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        for subj in meta:
            if subj.casefold() == "study": continue
            for group in meta[subj]["trials"]:                
                for trial in meta[subj]["trials"][group]:
                    trialoutpath = meta[subj]["trials"][group][trial]["outpath"]
                    if not os.path.exists(meta[subj]["outpath"]):
                        os.makedirs(meta[subj]["outpath"])
                    if not os.path.exists(os.path.join(meta[subj]["outpath"], group)):
                        os.makedirs(os.path.join(meta[subj]["outpath"], group))
                    if not os.path.exists(trialoutpath):
                        os.makedirs(trialoutpath)
                    shutil.copy(os.path.join(meta[subj]["trials"][group][trial]["inpath"], meta[subj]["trials"][group][trial]["c3dfile"]), trialoutpath)
                                        
    # save the metadata dict
    with open(os.path.join(outpath, user.project + ".pkl"), "wb") as fid:
        pk.dump(meta, fid)
    
    return meta


