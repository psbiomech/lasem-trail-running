# -*- coding: utf-8 -*-
"""
Write OpenSim data files and setup files

@author: Prasanna Sritharan
"""

import pandas as pd
import numpy as np
import os


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
    Write .mot file for ground forces (external loads).
'''
def write_ground_forces_mot_file(osimkey):

    # output dataframe info
    ns = len(osimkey.forces["time"])
    nc = 19    
    t0 = osimkey.forces["time"][0]
    tf = osimkey.forces["time"][-1]

    # write headers
    fname = osimkey.name + "_grf.mot"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("%s\n" % fname)
        f.write("nRows=%d\n" % ns)
        f.write("nColumns=%s\n" % nc)
        f.write("\n")
        f.write("name %s\n" % fname)
        f.write("datacolumns %d\n" % nc)
        f.write("datarows %d\n" % ns)
        f.write("range %f %f\n" % (t0, tf))
        f.write("endheader\n")

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



'''
write_marker_trajctory_trc_file(osimkey):
    Write .trc file for marker trajectories.
'''
def write_marker_trajctory_trc_file(osimkey):
    
    # output dataframe info
    ns = len(osimkey.markers["time"])
    nm = len(osimkey.markers) - 3
    rate = osimkey.markers["rate"]
    units = osimkey.markers["units"]

    # remove non-marker dict keys
    markernames0 = list(osimkey.markers.keys())
    markernames0.remove("rate")
    markernames0.remove("units")
    markernames0.remove("offset")
    markernames0.remove("frames")
    markernames0.remove("time")

    # build marker headers
    markernames = ""
    markernames = markernames.join(["Frame#\t","Time\t"] + list(map(lambda x: x + "\t\t\t", markernames0)))
    dirnums = list(range(1,len(markernames0) + 1))
    dirnames = ""
    dirnames = dirnames.join(["\t"] + list(map(lambda n: "\tX%d\tY%d\tZ%d" % (n, n, n), dirnums)))

    # write headers
    fname = osimkey.name + "_markers.trc"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("PathFileType\t4\t(X/Y/Z)\t%s\n" % fname)
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        f.write("%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (rate, rate, ns, nm, units, rate, 1, ns))
        f.write("%s\n" % markernames)
        f.write("%s\n" % dirnames)
        


    