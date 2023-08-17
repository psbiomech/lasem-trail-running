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

dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    traildb = pk.load(fid)


# %% LOAD AN OSIMKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL001_HFD_LEFT02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\hfd\hfd\TRAIL001\BASELINE\TRAIL001_HFD_LEFT02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_osimkey.pkl")
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
fprefix = "TRAIL001_EP02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stance\TRAIL001\BASELINE\TRAIL001_EP02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_c3dkey.pkl")
with open(pkfile, "rb") as fid: 
    c3dkey1 = pk.load(fid)
    
    


# %% LOAD AN OSIMRESULTSKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL006_HFD_LEFT03"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\hfd\hfd\TRAIL006\BASELINE\TRAIL006_HFD_LEFT03"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_opensim_results.pkl")
with open(pkfile, "rb") as fid: 
    osimresultskey1 = pk.load(fid)
    
    

# %% CALCULATE JOINT ANGULAR IMPULSE

import analyses as an

impl = an.calculate_joint_angular_impulse(osimkey1, user)