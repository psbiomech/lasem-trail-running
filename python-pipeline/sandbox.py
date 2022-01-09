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

# metadata
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
    
    
# %% RUN OPENSIM PIPELINE

import opensimpipeline as osp
import usersettings as uset
import pickle as pk
import os

# load user settings
user = uset.TRAILSettings()

# metadata
fpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data\outputDatabase"
fname = "TRAIL.pkl"
with open(os.path.join(fpath, fname),"rb") as fid: 
    traildb = pk.load(fid)

# run OpenSim pipeline
analyses = ["scale","ik","id"]
osp.opensim_pipeline(traildb, user, analyses)


    
# %% RUN SCALE TOOL

# assumption: OsimKey for static trial loaded

import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_Static_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey0 = pk.load(fid)

# run scale tool
osp.run_opensim_scale(osimkey0, user)



# %% RUN IK TOOL


import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey1 = pk.load(fid)

# run scale tool
osp.run_opensim_ik(osimkey1, user)



# %% RUN ID TOOL


import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey1 = pk.load(fid)

# run scale tool
osp.run_opensim_id(osimkey1, user)


# %% GET RESULTS


import opensimresults as osr
import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile,"rb") as fid: 
    osimkey1 = pk.load(fid)

# get results of dynamic trial using OsimKey
analyses = ["scale","ik","id"]
osimresult1 = osr.OsimResultsKey(osimkey1, analyses)

