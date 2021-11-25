# -*- coding: utf-8 -*-
"""
Write OpenSim data files and setup files

@author: Prasanna Sritharan
"""


import c3dextract as c3dex
import pandas as pd
import numpy as np


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
write_ground_forces_mot_file(osimkey):
    Write .mot file for ground forces (external loads)
'''
def write_ground_forces_mot_file(osimkey):

    # output dataframe info
    ns = len(osimkey.forces["time"])
    nc = 19    
    t0 = osimkey.forces["time"][0]
    tf = osimkey.forces["time"][-1]

    # write headers
    fname = osimkey.name + "_grf.mot"
    f = open(fname,"w")
    f.write("%s\n" % fname)
    f.write("nRows=%d\n" % ns)
    f.write("nColumns=%s\n" % nc)
    f.write("\n")
    f.write("name %s\n" % fname)
    f.write("datacolumns %d\n" % nc)
    f.write("datarows %d\n" % ns)
    f.write("range %f %f\n" % (t0, tf))
    f.write("endheader\n")
    f.close()

    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.forces["time"]
    datamat[:,1:4] = osimkey.forces["data"]["right"]["F"]
    datamat[:,4:7] = osimkey.forces["data"]["right"]["cop"]
    datamat[:,7:10] = osimkey.forces["data"]["right"]["T"]
    datamat[:,10:13] = osimkey.forces["data"]["left"]["F"]
    datamat[:,13:16] = osimkey.forces["data"]["left"]["cop"]
    datamat[:,16:19] = osimkey.forces["data"]["left"]["T"]    
        
    # convert to dataframe
    headers = ["time", "grf_right_vx", "grf_right_vy", "grf_right_vz", "grf_right_px", "grf_right_py", "grf_right_pz", "grf_left_vx", "grf_left_vy", "grf_left_vz", "grf_left_px", "grf_left_py", "grf_left_pz", "grf_right_tx", "grf_right_ty", "grf_right_tz", "grf_left_tx", "grf_left_ty", "grf_left_tz"]
    data = pd.DataFrame(datamat, columns=headers)
        
    # write table
    data.to_csv(fname, mode="a", sep="\t", header=True, index=False)
    
    return data
    