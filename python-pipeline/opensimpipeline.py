# -*- coding: utf-8 -*-
"""
Run OpenSim pipeline, write input data files

@author: Prasanna Sritharan
"""

import opensim
import pandas as pd
import numpy as np
import pickle as pk
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
run_opensim_scale(osimkey,user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run 
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_scale(osimkey, user):
    
    # reference setup file
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupscale
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    model = osimkey.subject
    trial = osimkey.trial

    print('Creating scaled model: %s\n' % model);
    
    # create an ScaleTool from a generic setup file
    tool = opensim.ScaleTool(os.path.join(refsetuppath, refsetupfile))
    tool.setPathToSubject("")
    
    # set subject mass
    tool.setSubjectMass(0)




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
    fname = osimkey.trial + "_grf.mot"
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
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=True, index=False)
    
    return data



'''
write_marker_trajctory_trc_file(osimkey):
    Write .trc file for marker trajectories.
'''
def write_marker_trajctory_trc_file(osimkey):
    
    # output dataframe info
    ns = len(osimkey.markers["time"])
    nm = len(osimkey.markers) - 5
    nc = 2 + (nm * 3)
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
    fname = osimkey.trial + "_markers.trc"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("PathFileType\t4\t(X/Y/Z)\t%s\n" % fname)
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        f.write("%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (rate, rate, ns, nm, "mm", rate, 1, ns))
        f.write("%s\n" % markernames)
        f.write("%s\n" % dirnames)
        
    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.markers["frames"]
    datamat[:,1] = osimkey.markers["time"]
    n = 2
    for mkr in markernames0:
        mkrdata = osimkey.markers[mkr]
        if units.casefold() == "m": mkrdata = mkrdata * 1000
        datamat[:,n:n+3] = mkrdata
        n = n + 3
    
    # convert to dataframe
    data = pd.DataFrame(datamat)
    data[0] = data[0].astype(int)
    
    # write table, no headers
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=False, index=False)
    
    return data