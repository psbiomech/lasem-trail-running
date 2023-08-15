# -*- coding: utf-8 -*-
"""
Copy C3D files from P:/ network drive to local TRAIL/inputDatabase.

Note: must be logged into network drive.

@author: Prasanna Sritharan
"""

import os
import glob
import re
import shutil


pfolder = r"P:\PROJECT - OPV - TRAIL\TRAIL\Gait Lab and Biodex\Data-Biomech\4a Checked and ready_running and HFD"
cfolder = r"C:\Users\Owner\Documents\data\TRAIL\inputdatabase"


# get list of C3D files
inpath = os.path.join(pfolder, "**", "*.c3d")
folderlist = glob.glob(inpath, recursive=True)

# file name expression
#fstr = "TRAIL(?:_?)(\d+)_(EP|FAST|STATIC)_?(\d+)"
fstr = "TRAIL(?:_?)(\d+)_(HFD_LEFT|HFD_RIGHT|STATIC)_?(\d+)"
fregx = re.compile(fstr, re.IGNORECASE)

# get the file name and parse with regex
if not(os.path.isdir(cfolder)): os.makedirs(cfolder)
for f in folderlist:
    
    # get the subject and trial code
    fname0 = os.path.splitext(os.path.split(f)[1])
    fcodes = fregx.match(fname0[0].upper())
    if not fcodes: continue
    
    # create new folders
    subjdir = os.path.join(cfolder, "TRAIL" + fcodes[1])
    if not(os.path.isdir(subjdir)): os.makedirs(subjdir)
    phasedir = os.path.join(subjdir, "BASELINE")
    if not(os.path.isdir(phasedir)): os.makedirs(phasedir)
    
    # copy C3D file to new folder
    fname1 = "TRAIL" + fcodes[1] + "_" + fcodes[2].upper() + fcodes[3] + ".c3d"
    shutil.copy(f, os.path.join(phasedir, fname1))
    







