# -*- coding: utf-8 -*-
"""
Main: extract TRAIL test data

@author: Prasanna Sritharan
"""

#import numpy as np
import c3dextract as c3dex
import labsetup as labs


# set the lab
lasem = labs.LabKeyLasemTrail()

# return C3DKey and TrialKey
c3dfile = r"C:\Users\Owner\Documents\data\TRAIL Test Data\Events\TRAIL_071\Baseline\TRAIL_071_EP_01.c3d"
c3dkey, trialkey, osimkey = c3dex.c3d_extract(c3dfile, lasem, 2, 15)

