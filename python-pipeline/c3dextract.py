# -*- coding: utf-8 -*-
"""
LASEM C3D file data extract
@author: Prasanna Sritharan

"""

import numpy as np
import pyc3dserver as c3d



'''
---------------------------------
------------ CLASSES ------------
---------------------------------
'''

'''
C3DKey:
    C3D data storage class containing all C3D file data
'''
class C3DKey():
    def __init__(self, fname, fmeta, fforces, fmarkers):
        self.name = fname
        self.meta = fmeta
        self.forces = fforces
        self.markers = fmarkers



'''
TrialKey:
    Trial data storage class containing all C3D data relevant to trial in the
    laboratory and force plate frames
'''
class TrialKey():
    def __init__(self,lab,c3dkey):
        self.markers = self.set_markers(c3dkey)
        self.forces = self.set_forces(lab,c3dkey)
        

    def set_markers(self,c3dkey):
        
        # initialise dict
        markers = {}
        
        # marker rate and units
        markers["rate"] = c3dkey.markers["RATE"]    
        markers["units"] = c3dkey.markers["UNITS"]
        
        # marker scale factor
        if markers["units"] == "mm":
            markers["scale"] = 0.001
        elif markers["units"] == "m":
            markers["scale"] = 1
        else:
            markers["scale"] = 1

        # marker data
        markers["labels"] = c3dkey.markers["LABELS"]
        markers["data"] = c3dkey.markers["DATA"]["POS"]
        markers["time"] = c3dkey.markers["TIME"]
        markers["time0"] = c3dkey.markers["TIME"] - c3dkey.markers["TIME"][0]    
        markers["frames"] = c3dkey.markers["FRAME"]
        markers["frames0"] = c3dkey.markers["FRAME"] - c3dkey.markers["FRAME"][0]    
        
        return markers
           
    def set_forces(self,lab,c3dkey):
        
        # initialise dict
        forces = {}        
 
        # get forces for only used force plate
        for f in lab.fpused:
            
            # create force plate channel names (lazy formatting)
            channels = [s + str(f) for s in lab.fpprefixes]
            channels = [i + j for i, j in zip(channels, lab.fpsuffixes)]
     
            # force plate labels: find first index, assume the others are in
            # the correct order
            indx0 = np.where(c3dkey.forces["LABELS"]==channels[0])[0][0]
            indx = list(range(indx0,indx0+lab.nfpchannels))
            forces["labels"] = np.ndarray.tolist(np.array(c3dkey.forces["LABELS"])[indx])
     
            # force rate and units
            forces["rate"] = c3dkey.forces["RATE"]    
            forces["units"] = np.ndarray.tolist(np.array(c3dkey.forces["UNITS"])[indx])
            
            # force scale factor
            forces["scale"] = []
            for x in range(lab.nfpchannels):
                if forces["units"][x] == "Nmm":
                    forces["scale"].append(0.001)
                elif forces["units"][x] == "Nm":
                    forces["scale"].append(1)
                elif forces["units"][x] == "N":
                    forces["scale"].append(1)
                else:
                    forces["scale"].append(1)
            
            # force plate time and frames
            forces["time"] = c3dkey.forces["TIME"]
            forces["time0"] = c3dkey.forces["TIME"] - c3dkey.forces["TIME"][0]    
            forces["frames"] = c3dkey.forces["FRAME"]
            forces["frames0"] = c3dkey.forces["FRAME"] - c3dkey.forces["FRAME"][0] 
        
            # force plate data
            forces["data"] = {}
            for ch in channels:
                forces["data"][ch] = c3dkey.forces["DATA"][ch]
        
       
        return forces
    


'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
c3d_extract(f_path):
    Extracts the motion data from the C3D file to arrays, and returns a dict
    containing all the relevant file metadata, force data and marker data.
'''
def c3d_extract(f_path,lab):
    
    # load C3D file
    itf = c3d.c3dserver()
    c3d.open_c3d(itf, f_path)       
    
    # get all file metadata, and all force plate and video C3D data
    fmeta = c3d.get_dict_groups(itf)
    fforces = c3d.get_dict_forces(itf, frame=True, time=True)
    fmarkers = c3d.get_dict_markers(itf, frame=True, time=True)
    
    # subject name
    fname= fmeta["SUBJECTS"]["NAMES"][0]
    if not fname:
        fname = "NoName"
    
    # C3D key with all data from C3D file
    c3dkey = C3DKey(fname, fmeta, fforces, fmarkers)

    # trial data only from C3D key
    trialkey = TrialKey(lab,c3dkey)
    
    return c3dkey, trialkey
    

'''
calc_cop(meta,forces):
    calculate centre-of-pressure (CoP)
'''
#def calc_cop(meta,forces):
