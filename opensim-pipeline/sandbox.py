# -*- coding: utf-8 -*-
"""
Sandbox

@author: Prasanna Sritharan
"""


# %% LOAD USER SETTINGS

import usersettings as uset

user = uset.FORCESettings_SDP()



# %% LOAD META DATA

import pickle as pk
import os

dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    forcedb = pk.load(fid)


# %% LOAD AN OSIMKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL022_STATIC02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\TRAIL001\BASELINE\TRAIL001_EP02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_osimkey.pkl")
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)
 
    
# %% LOAD A TRIALKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL022_EP03"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\TRAIL022\BASELINE\TRAIL022_EP03"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_trialkey.pkl")
with open(pkfile, "rb") as fid: 
    trialkey1 = pk.load(fid)    
    
# %% LOAD A C3DKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL022_STATIC02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\TRAIL006\BASELINE\TRAIL006_STATIC02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_c3dkey.pkl")
with open(pkfile, "rb") as fid: 
    c3dkey1 = pk.load(fid)
    
    


# %% LOAD AN OSIMRESULTSKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "TRAIL001_FAST02"
fpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\TRAIL001\BASELINE\TRAIL001_FAST02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_opensim_results.pkl")
with open(pkfile, "rb") as fid: 
    osimresultskey1 = pk.load(fid)
    
    

# %% CALCULATE JOINT ANGULAR IMPULSE

import analyses as an

impl = an.calculate_joint_angular_impulse(osimkey1, user)