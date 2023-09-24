# -*- coding: utf-8 -*-
"""
Sandbox

@author: Prasanna Sritharan
"""


# %% LOAD USER SETTINGS

import usersettings as uset

user = uset.TRAILSettings_RUN()



# %% LOAD META DATA

import pickle as pk
import os

dbfilepath = os.path.join(user.rootpath, user.outfolder,"run", "run_stridecycle", user.metadatafile)
with open(dbfilepath, "rb") as fid:
    traildb = pk.load(fid)


# %% LOAD AN OSIMKEY

import pickle as pk
import os

# file path and name prefix
subj = "TRAIL225"
trial = "EP08"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stridecycle"

# OsimResultsKey
pkfile = os.path.join(fpath, subj, "BASELINE", subj + "_" + trial, subj + "_" + trial + "_osimkey.pkl")
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)
 
    
# %% LOAD A TRIALKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL001_EP02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stance\TRAIL001\BASELINE\TRAIL001_EP02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_trialkey.pkl")
with open(pkfile, "rb") as fid: 
    trialkey1 = pk.load(fid)    
    
# %% LOAD A C3DKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL123_EP02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stridecycle\TRAIL123\BASELINE\TRAIL123_EP02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_c3dkey.pkl")
with open(pkfile, "rb") as fid: 
    c3dkey1 = pk.load(fid)
    
    

# %% LOAD AN OSIMRESULTSKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL001_EP02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stance\TRAIL001\BASELINE\TRAIL001_EP02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_opensim_results.pkl")
with open(pkfile, "rb") as fid: 
    osimresultskey1 = pk.load(fid)
    
    

# %% TEST ANALYSES

import usersettings as uset
import analyses as an

user = uset.TRAILSettings_RUN()

jap = an.calculate_joint_angular_power(osimresultskey1, user)
jaw = an.calculate_joint_angular_work(osimresultskey1, user)




