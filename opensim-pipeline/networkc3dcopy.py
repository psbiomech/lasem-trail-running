# -*- coding: utf-8 -*-
"""
Copy C3D files from network drive to inputDatabase

@author: Prasanna Sritharan
"""

import os
import glob
import re
import shutil



pfolder = r"P:\PROJECT - SHE - ACLKNEE\TRAIL\Gait Lab and Biodex\Data-Biomech\3 FloorEvents_Added"
cfolder = r"C:\Users\Owner\Documents\data\TRAIL\inputDatabase"



# get list of C3D files
inpath = os.path.join(pfolder, "**", "*.c3d")
folderlist = glob.glob(inpath, recursive=True)

# file name expression
fstr = "TRAIL_(\d+)_(EP|FAST|STATIC|Static)_?(\d+)"
fregx = re.compile(fstr)

# get the file name and parse with regex
if not(os.path.isdir(cfolder)): os.makedirs(cfolder)
for f in folderlist:
    
    # get the subject and trial code
    fname0 = os.path.splitext(os.path.split(f)[1])
    fcodes = fregx.match(fname0[0])
    if not fcodes: continue
    
    # create new folders
    subjdir = os.path.join(cfolder, "TRAIL" + fcodes[1])
    if not(os.path.isdir(subjdir)): os.makedirs(subjdir)
    phasedir = os.path.join(subjdir, "BASELINE")
    if not(os.path.isdir(phasedir)): os.makedirs(phasedir)
    
    # copy C3D file to new folder
    fname1 = "TRAIL" + fcodes[1] + "_" + fcodes[2].upper() + fcodes[3] + ".c3d"
    shutil.copy(f, os.path.join(phasedir, fname1))
    







