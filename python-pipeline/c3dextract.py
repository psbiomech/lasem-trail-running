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
    def __init__(self, lab, c3dkey, xdir):
        
        self.trial_name = str(c3dkey.name)
        self.lab_name = lab.lab_name
        
        self.__set_force_plates(lab, c3dkey, xdir) 
        self.__set_forces(lab,c3dkey)
        self.__set_markers(c3dkey)

        
    def __set_markers(self, c3dkey):
        
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
        
        self.markers = markers
        
        return None
           
    def __set_forces(self, lab, c3dkey):
        
        # initialise dict
        forces = {}        
 
        # get forces for only used force plates
        for f in lab.fp_used:
            
            # dict field
            dict_name = lab.fp_dict_name_prefix + str(f)
            forces[dict_name] = {}
            
            # create force plate channel names (lazy formatting)
            channels = [s + str(f) for s in lab.fp_channel_prefixes]
            channels = [i + j for i, j in zip(channels, lab.fp_channel_suffixes)]
     
            # force plate labels: find first index, assume the others are in
            # the correct order
            indx0 = np.where(c3dkey.forces["LABELS"]==channels[0])[0][0]
            indx = list(range(indx0,indx0+lab.n_fp_channels))
            forces[dict_name]["labels"] = np.ndarray.tolist(np.array(c3dkey.forces["LABELS"])[indx])
     
            # force rate and units
            forces[dict_name]["rate"] = c3dkey.forces["RATE"]    
            forces[dict_name]["units"] = np.ndarray.tolist(np.array(c3dkey.forces["UNITS"])[indx])
            
            # force scale factor
            forces[dict_name]["scale"] = []
            for x in range(lab.n_fp_channels):
                if forces[dict_name]["units"][x] == "Nmm":
                    forces[dict_name]["scale"].append(0.001)
                elif forces[dict_name]["units"][x] == "Nm":
                    forces[dict_name]["scale"].append(1)
                elif forces[dict_name]["units"][x] == "N":
                    forces[dict_name]["scale"].append(1)
                else:
                    forces[dict_name]["scale"].append(1)
            
            # force plate time and frames
            forces[dict_name]["time"] = c3dkey.forces["TIME"]
            forces[dict_name]["time0"] = c3dkey.forces["TIME"] - c3dkey.forces["TIME"][0]    
            forces[dict_name]["frames"] = c3dkey.forces["FRAME"]
            forces[dict_name]["frames0"] = c3dkey.forces["FRAME"] - c3dkey.forces["FRAME"][0] 
        
            # force plate data
            forces[dict_name]["data"] = {}
            for ch in channels:
                forces[dict_name]["data"][ch] = c3dkey.forces["DATA"][ch]
        
       
        self.forces = forces
        
        return None
    
    def __set_force_plates(self, lab, c3dkey, xdir):
        
        # initialise dict
        force_plates = {}
        
        # get force plate info for only used force plates
        for f in lab.fp_used:
            
            # dict field
            dict_name = lab.fp_dict_name_prefix + str(f)
            force_plates[dict_name] = {}
            
            # coordinate transforms
            force_plates[dict_name]["transforms"] = {}
            force_plates[dict_name]["transforms"]["lab_to_opensim"] = lab.transform_mat_lab_to_opensim[[1, -1, 2, -2, 3, -3].index(xdir)]
            force_plates[dict_name]["transforms"]["fp_to_lab"] = lab.transform_mat_fp_to_lab
        
            # offsets
            force_plates[dict_name]["offsets"] = {}
            force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"] = -1*c3dkey.meta["FORCE_PLATFORM"]["ORIGIN"][f-1]
            force_plates[dict_name]["offsets"]["lab_to_fp_centre_vicon"] = find_fp_centre_from_lab_origin(c3dkey.meta["FORCE_PLATFORM"]["CORNERS"][f-1])
            force_plates[dict_name]["offsets"]["lab_to_fp_origin_vicon"] = change_coordinates(force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"], force_plates[dict_name]["transforms"]["fp_to_lab"], force_plates[dict_name]["offsets"]["lab_to_fp_centre_vicon"])
        
        self.force_plates = force_plates
        
        return None



'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
c3d_extract(f_path, lab, xdir):
    Extract the motion data from the C3D file to arrays, and returns a dict
    containing all the relevant file metadata, force data and marker data.
'''
def c3d_extract(f_path,lab,xdir):
    
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
    trialkey = TrialKey(lab, c3dkey, xdir)
    
    return c3dkey, trialkey
    



'''
find_fp_centre_from_lab_origin(corners):
    Find force plate centre in Vicon coordinates
'''
def find_fp_centre_from_lab_origin(corners):
    
    # lab origin ---> corner 1
    lo_c1 = corners[0]
    
    # lab origin ---> corner 3
    lo_c3 = corners[2]
    
    # corner 1 ---> corner 3
    c1_c3 = lo_c3 - lo_c1
    
    # corner 1 ---> centre
    c1_ct = c1_c3 / 2
    
    # lab origin ---> centre
    return lo_c1 + c1_ct




'''
change_coordinates(oldvec, rotmat, originvec):
    Perform coordinate transformation
'''
def change_coordinates(oldvec, rotmat, originvec):
    return originvec + np.dot(rotmat, oldvec)