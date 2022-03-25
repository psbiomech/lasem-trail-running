# -*- coding: utf-8 -*-
"""
Sandbox

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
with open(pkfile, "rb") as fid: 
    traildb = pk.load(fid)



# %% LOAD A STATIC TRIAL

import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_Static_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# C3DKey
pkfile = prefix + "_c3dkey.pkl"
with open(pkfile, "rb") as fid: 
    c3dkey0 = pk.load(fid)

# TrialKey
pkfile = prefix + "_trialkey.pkl"
with open(pkfile, "rb") as fid: 
    trialkey0 = pk.load(fid)
    
# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey0 = pk.load(fid)
    


# %% LOAD A DYNAMIC TRIAL

import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_FAST_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# C3DKey
pkfile = prefix + "_c3dkey.pkl"
with open(pkfile, "rb") as fid: 
    c3dkey1 = pk.load(fid)

# TrialKey
pkfile = prefix + "_trialkey.pkl"
with open(pkfile, "rb") as fid: 
    trialkey1 = pk.load(fid)
    
# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
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
with open(os.path.join(fpath, fname), "rb") as fid: 
    traildb = pk.load(fid)

# run OpenSim pipeline
analyses = ["scale", "ik", "id"]
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
with open(pkfile, "rb") as fid: 
    osimkey0 = pk.load(fid)

# run tool
osp.run_opensim_scale(osimkey0, user)



# %% RUN IK TOOL


import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_FAST_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# run tool
osp.run_opensim_ik(osimkey1, user)



# %% RUN ID TOOL


import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_FAST_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# run tool
osp.run_opensim_id(osimkey1, user)


# %% RUN SO TOOL


import pickle as pk
import opensimpipeline as osp
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# file path and name prefix
fprefix = "TRAIL_071_FAST_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# run tool
osp.run_opensim_so(osimkey1, user)


# %% RUN RRA TOOL


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
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# run tool
osp.run_opensim_rra(osimkey1, user)


# %% RUN CMC TOOL


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
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# run tool
osp.run_opensim_cmc(osimkey1, user)


# %% CREATE OSIMRESULTSKEY FROM OSIMKEY AND OPENSIM RESULTS


import opensimresults as osr
import pickle as pk

# file path and name prefix
#fprefix = "TRAIL_071_FAST_02"
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)

# get results of dynamic trial using OsimKey
analyses = ["scale", "ik", "id", "so"]
osimresult1 = osr.OsimResultsKey(osimkey1, analyses, 101)


# %% RESAMPLE TEST


import numpy as np
import opensimresults as osr

test0  = np.array([[1, 2, 3, 4, 5],[2, 4, 6, 8, 10]]).transpose()
test1 = osr.resample1d(test0, 11)


# %% OSIMRESULTSKEY BATCH PROCESS TEST

import pickle as pk
import opensimresults as osr
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# metadata
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\"
fname = "TRAIL.pkl"
pkfile = fpath + fname
with open(pkfile, "rb") as fid: 
    traildb = pk.load(fid)

# batch process
analyses = ["so"]
osr.opensim_results_batch_process(traildb, analyses, 101)


# %% LOAD AN OSIMRESULTSKEY

import pickle as pk

# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimResultsKey
pkfile = prefix + "_opensim_results.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)



# %% COLLATE AND EXPORT TEST


import pickle as pk
import opensimresults as osr
import usersettings as uset

# load user settings
user = uset.TRAILSettings()

# metadata
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\"
fname = "TRAIL.pkl"
pkfile = fpath + fname
with open(pkfile,"rb") as fid: 
    traildb = pk.load(fid)

# collate and export
analyses = ["ik", "id"]
csvdata = osr.export_opensim_results(traildb, user, analyses)


# %% FILTERING DATA

import pickle as pk
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


# file path and name prefix
fprefix = "TRAIL_071_EP_01"
fpath = "C:\\Users\\Owner\\Documents\\data\\TRAIL Test Data\\outputDatabase\\TRAIL_071\\Baseline\\" + fprefix + "\\"
prefix = fpath + fprefix

# OsimKey
pkfile = prefix + "_osimkey.pkl"
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)
    
# data
grf0 = osimkey1.forces["data"]["left"]["F"][:, 1]

# filter design
samplerate = 1200.0
Wn = samplerate / 2
cutoff = 50
normalised_cutoff = cutoff / Wn
b, a = signal.butter(4, normalised_cutoff, "lowpass")

# filtered data
grf1 = signal.filtfilt(b, a, grf0)

# floor below threshold
threshold = 10
idxs = np.where(grf1 < threshold)
grf1[idxs] = 0

# plot
grfs = np.transpose([grf0, grf1])
plt.plot(grfs)
