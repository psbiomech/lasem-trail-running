# -*- coding: utf-8 -*-
"""
OpenSim pipeline sandbox

@author: Prasanna Sritharan
"""



# %% LOAD USER SETTINGS

import usersettings as uset

# load user settings
user = uset.TRAILSettings()





# %% LOAD METADATA DICT

import pickle as pk

# file path and name prefix
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\"
fname = "TRAIL.pkl"
pkfile = fpath + fname
with open(pkfile,"rb") as fid: 
    traildb = pk.load(fid)



# %% LOAD A STATIC TRIAL

import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_Static_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# C3DKey
pkfile = prefix + "_c3dkey.pkl"
with open(pkfile,"rb") as fid: 
    c3dkey0 = pk.load(fid)

# TrialKey
pkfile = prefix + "_trialkey.pkl"
with open(pkfile,"rb") as fid: 
    trialkey0 = pk.load(fid)
    
# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey0 = pk.load(fid)
    


# %% LOAD A DYNAMIC TRIAL

import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# C3DKey
pkfile = prefix + "_c3dkey.pkl"
with open(pkfile,"rb") as fid: 
    c3dkey1 = pk.load(fid)

# TrialKey
pkfile = prefix + "_trialkey.pkl"
with open(pkfile,"rb") as fid: 
    trialkey1 = pk.load(fid)
    
# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey1 = pk.load(fid)
    
    
