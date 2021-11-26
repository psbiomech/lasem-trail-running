# -*- coding: utf-8 -*-
"""
LASEM C3D file data extract

@author: Prasanna Sritharan
"""

import numpy as np
import pyc3dserver as c3d
import pickle as pk
import os



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

'''
C3DKey:
    C3D data storage class containing all C3D file data.
'''
class C3DKey():
    def __init__(self, fname, fmeta, fforces, fmarkers):
        self.name = fname
        self.meta = fmeta
        self.forces = fforces
        self.markers = fmarkers
        return None



'''
TrialKey:
    Trial data storage class containing all C3D data relevant to trial in the
    laboratory and force plate frames.
'''
class TrialKey():
    def __init__(self, lab, c3dkey, xdir):       
        self.trial_name = str(c3dkey.name)
        self.lab_name = lab.lab_name
        self.__set_events(c3dkey)
        self.__set_markers(lab, c3dkey, xdir)        
        self.__set_force_plates(lab, c3dkey, xdir) 
        self.__set_forces(lab, c3dkey)
        return None

    def __set_events(self, c3dkey):

        # initialise dict        
        events = {}
        
        # return empty dict if no events
        if c3dkey.meta["EVENT"]["USED"] == 0:
            self.events = events
            return None
        
        # all event labels
        foot = [f[0].upper() for f in c3dkey.meta["EVENT"]["CONTEXTS"]]
        events["labels"] = [foot[i] + "F" + f.split()[1][0] for i, f in enumerate(c3dkey.meta["EVENT"]["LABELS"])]
        
        # all event times
        events["time"] = c3dkey.meta["EVENT"]["TIMES"][:,1]
        events["time0"] = c3dkey.meta["EVENT"]["TIMES"][:,1] - (c3dkey.meta["TRIAL"]["ACTUAL_START_FIELD"][0] / c3dkey.meta["TRIAL"]["CAMERA_RATE"])
        
        
    
        
        # *********************
        # RUNNING ONLY (HARD-CODED TEMPORARILY)
        
        # temporarily assume running task only, so only process first full
        # stride cycle (IFS ---> CFO), will need to change this later to accept
        # any sequence to allow for processing different types of tasks
        
        # calculate the time window of interest (assume first FS is start of
        # the trial time window
        fsidx0 = np.where(np.char.find(events["labels"],"FS")>=0)[0][0]
        foidx1 = np.where(np.char.find(events["labels"],"FO")>=0)[0][1]
        events["window_time0"] = events["time0"][fsidx0:foidx1 + 1]
        events["window_labels"] = events["labels"][fsidx0:foidx1 + 1]
        
        # list the individual intervals
        events["window_intervals0"] = np.array([[t0, t1] for t0, t1 in zip(events["window_time0"][0:-1], events["window_time0"][1:])])
        
        # find force plate sequence for each interval (row) defined in the
        # array window_intervals0
        # note: sequences defined as per GaitExtract using a 2D array:
        #   rows: event intervals
        #   col1: right foot
        #   col2: left foot
        if events["window_labels"][0][0] == "R":
            events["fp_sequence"] = np.array([[4, 0], [0, 0], [0, 3]])
        else:
            events["fp_sequence"] = np.array([[0, 4], [0, 0], [3, 0]],)
        
        # *********************
        
        
        
        self.events = events
        
        return None
        
    def __set_markers(self, lab, c3dkey, xdir):
        
        # initialise dict
        markers = {}
        markers["labels"] = c3dkey.markers["LABELS"]
        
        # marker rate
        markers["rate"] = c3dkey.markers["RATE"]    
        orig_units = c3dkey.markers["UNITS"]
        
        # marker scale factor, convert to m
        markers["units"] = "m"
        if orig_units == "mm":
            markers["scale"] = 0.001
        elif orig_units == "m":
            markers["scale"] = 1
        else:
            markers["scale"] = 1

        # coordinate transforms
        markers["transforms"] = {}
        markers["transforms"]["lab_to_opensim"] = lab.transform_mat_lab_to_opensim[[1, -1, 2, -2, 3, -3].index(xdir)]

        # offset.offset_marker marker
        markers["offset_marker"] = lab.offset_marker
        
        # marker time and frames
        markers["time"] = c3dkey.markers["TIME"]
        markers["time0"] = c3dkey.markers["TIME"] - c3dkey.markers["TIME"][0]    
        markers["frames"] = c3dkey.markers["FRAME"]
        markers["frames0"] = c3dkey.markers["FRAME"] - c3dkey.markers["FRAME"][0]    

        # marker data, convert to m
        markers["data"] = {}
        for mkr in c3dkey.markers["DATA"]["POS"].keys():
            markers["data"][mkr] = c3dkey.markers["DATA"]["POS"][mkr] * markers["scale"]
        
        self.markers = markers
        
        return None
     
    def __set_force_plates(self, lab, c3dkey, xdir):
        
        # initialise dict
        force_plates = {}
        
        # get force plate info for only used force plates
        force_plates["fp_used"] = [] 
        force_plates["fp_used_str"] = []
        for f in lab.fp_used:
            
            # dict field
            dict_name = lab.fp_dict_name_prefix + str(f)
            force_plates["fp_used"].append(f)
            force_plates["fp_used_str"].append(dict_name)          
            
            # coordinate transforms
            force_plates[dict_name] = {}  
            force_plates[dict_name]["transforms"] = {}
            force_plates[dict_name]["transforms"]["lab_to_opensim"] = lab.transform_mat_lab_to_opensim[[1, -1, 2, -2, 3, -3].index(xdir)]
            force_plates[dict_name]["transforms"]["fp_to_lab"] = lab.transform_mat_fp_to_lab 
        
            # offsets
            offset_scale = self.markers["scale"]
            force_plates[dict_name]["offsets"] = {}
            force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"] = -1*c3dkey.meta["FORCE_PLATFORM"]["ORIGIN"][f-1] * offset_scale
            force_plates[dict_name]["offsets"]["lab_to_fp_centre_lab"] = find_fp_centre_from_lab_origin(c3dkey.meta["FORCE_PLATFORM"]["CORNERS"][f-1]) * offset_scale
            force_plates[dict_name]["offsets"]["lab_to_fp_origin_lab"] = change_coordinates(force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"], force_plates[dict_name]["transforms"]["fp_to_lab"], force_plates[dict_name]["offsets"]["lab_to_fp_centre_lab"])
        
        self.force_plates = force_plates
        
        return None
          
    def __set_forces(self, lab, c3dkey):
        
        # initialise dict
        forces = {}        
              
        # force plate time and frames
        forces["time"] = c3dkey.forces["TIME"]
        forces["time0"] = c3dkey.forces["TIME"] - c3dkey.forces["TIME"][0]    
        forces["frames"] = c3dkey.forces["FRAME"]
        forces["frames0"] = c3dkey.forces["FRAME"] - c3dkey.forces["FRAME"][0]    
 
        # get forces for only used force plates
        ns = len(forces["time"])
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
            old_units = np.ndarray.tolist(np.array(c3dkey.forces["UNITS"])[indx])
            
            # force scale factor, convert to N and Nm
            forces[dict_name]["scale"] = []
            forces[dict_name]["units"] = []
            for x in range(lab.n_fp_channels):
                if old_units[x] == "Nmm":
                    forces[dict_name]["scale"].append(0.001)
                    forces[dict_name]["units"].append("Nm")
                elif old_units[x] == "Nm":
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("Nm")
                elif old_units[x] == "N":
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("N")
                else:
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("N")
        
            # force plate data, convert to N and Nm
            forces[dict_name]["data"] = {}
            F = np.zeros([ns,3])
            M = np.zeros([ns,3])
            for c in range(3):
                F[:,c] = c3dkey.forces["DATA"][channels[c]] * forces[dict_name]["scale"][c]
                M[:,c] = c3dkey.forces["DATA"][channels[c+3]] * forces[dict_name]["scale"][c+3]                
            forces[dict_name]["data"]["F"] = F
            forces[dict_name]["data"]["M"] = M
                
            # calculate centre-of-pressure
            vc2o = self.force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"]                        
            cop = calculate_centre_of_pressure_fp(ns, vc2o, F, M)
            forces[dict_name]["data"]["cop"] = cop
                
            # calculate free moments
            T = calculate_vertical_free_moment(ns, vc2o, F, M, cop)
            forces[dict_name]["data"]["T"] = T
            
            # convert to lab (Vicon) coordinates
            originvec0 = [0, 0, 0]
            rotmat = self.force_plates[dict_name]["transforms"]["fp_to_lab"] 
            originvec = self.force_plates[dict_name]["offsets"]["lab_to_fp_origin_lab"]                    
            F_lab = np.zeros([ns,3])
            M_lab = np.zeros([ns,3])
            cop_lab = np.zeros([ns,3])
            T_lab = np.zeros([ns,3])
            for n in range(ns):
                F_lab[n,:] = change_coordinates(F[n,:], rotmat, originvec0)
                M_lab[n,:] = change_coordinates(M[n,:], rotmat, originvec0)
                cop_lab[n,:] = change_coordinates(cop[n,:], rotmat, originvec)
                T_lab[n,:] = change_coordinates(T[n,:], rotmat, originvec0)
            forces[dict_name]["data"]["F_lab"] = F_lab
            forces[dict_name]["data"]["M_lab"] = M_lab
            forces[dict_name]["data"]["cop_lab"] = cop_lab
            forces[dict_name]["data"]["T_lab"] = T_lab
            
        self.forces = forces
        
        return None




'''
OpenSimKey:
    OpenSim data storage class for all data, model names, and other data for 
    processing through OpenSim.
'''
class OpenSimKey():
    def __init__(self, trialkey, ref_model, c3dpath, threshold):
        self.name = trialkey.trial_name
        self.lab = trialkey.lab_name
        self.model = trialkey.trial_name + ".osim"
        self.outpath = c3dpath
        self.ref_model = ref_model
        self.__set_markers(trialkey) 
        self.__set_forces(trialkey, threshold)      
        return None

    def __set_markers(self, trialkey):
        
        # initialise dict
        markers = {}
        
        # time and frame vectors
        markers["time"] = trialkey.markers["time0"]
        markers["frames"] = np.arange(1,len(trialkey.markers["time0"]) + 1)             
        
        # get markers
        data = {}
        for mkr in trialkey.markers["data"].keys():
            
            # current marker data
            data_lab = trialkey.markers["data"][mkr]
            
            # convert data from lab coordinates to OpenSim coordinates
            ns = len(trialkey.markers["time0"])
            originvec0 = [0, 0, 0]
            rotmat = trialkey.markers["transforms"]["lab_to_opensim"]
            data[mkr] = np.zeros([ns,3])
            for n in range(ns): data[mkr][n,:] = change_coordinates(data_lab[n,:], rotmat, originvec0)
                  
        # marker offset
        markers["offset"] = [0., 0., 0.]
        if trialkey.markers["offset_marker"]: markers["offset"] = [data[trialkey.markers["offset_marker"]][0, 0], 0., data[trialkey.markers["offset_marker"]][0, 2]]
                       
        # offset X and Z trajectories (OpenSim coordinates system) of markers
        data_offset = data
        for mkr in data.keys():
            for n in range(ns): data_offset[mkr][n,:] = data[mkr][n,:] - markers["offset"]
            markers[mkr] = data_offset[mkr]
            
        self.markers = markers
        
        return None
                   
    def __set_forces(self, trialkey, threshold):
 
        # initiliase temporary output arrays
        data = {}
        leg = ["right","left"]
        ns = len(trialkey.forces["time0"])
        for h, g in enumerate(leg):
            data[g] = {}
            data[g]["F"] = np.zeros([ns,3])
            data[g]["cop"] = np.zeros([ns,3])
            data[g]["T"] = np.zeros([ns,3])
    
        # initialise dict
        forces = {}        
      
        # time and frame vectors
        forces["time"] = trialkey.forces["time0"]
        forces["frames"] = np.arange(1,len(trialkey.forces["time0"]) + 1)

        # get forces for only used force plates
        forces["data"] = {}
        for i, fp in enumerate(trialkey.force_plates["fp_used"]):
            
            # dict field
            dict_name = trialkey.force_plates["fp_used_str"][i]

            # force data in lab (Vicon) coordinates
            F_lab = trialkey.forces[dict_name]["data"]["F_lab"]
            cop_lab = trialkey.forces[dict_name]["data"]["cop_lab"]
            T_lab = trialkey.forces[dict_name]["data"]["T_lab"]

            # convert data from lab coordinates to OpenSim coordinates
            originvec0 = [0, 0, 0]
            rotmat = trialkey.force_plates[dict_name]["transforms"]["lab_to_opensim"]
            F = np.zeros([ns,3])
            cop = np.zeros([ns,3])
            T = np.zeros([ns,3])
            for n in range(ns):
                F[n,:] = change_coordinates(F_lab[n,:], rotmat, originvec0)
                cop[n,:] = change_coordinates(cop_lab[n,:], rotmat, originvec0)
                T[n,:] = change_coordinates(T_lab[n,:], rotmat, originvec0)            
            
            # offset the CoP using offset marker
            offset = [0., 0., 0.]
            if trialkey.markers["offset_marker"]: offset = self.markers["offset"]
            for n in range(ns): cop[n,:] = cop[n,:] - offset

            # for each force plate, add force plate data for any active
            # intervals to the output array for the relevant foot
            for h, g in enumerate(leg):
                for n, m in enumerate(trialkey.events["fp_sequence"][:,h]):
                    if m == fp:                        
                        idx0 = np.where(forces["time"] >= trialkey.events["window_intervals0"][n,0])[0][0]
                        idx1 = np.where(forces["time"] <= trialkey.events["window_intervals0"][n,1])[0][-1]
                        data[g]["F"][idx0:idx1+1,:] = F[idx0:idx1+1,:]
                        data[g]["cop"][idx0:idx1+1,:] = cop[idx0:idx1+1,:]
                        data[g]["T"][idx0:idx1+1,:] = T[idx0:idx1+1,:]
            
        # store force data
        forces["data"] = data
        
        self.forces = forces
        
        return None




'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
c3d_batch_extract(meta, lab, xdir, threshold, ref_model):
    Batch processing for c3d_extract().
'''
def c3d_batch_extract(meta, lab, xdir, threshold, ref_model):
    for subj in meta:
        for group in meta[subj]["trials"]:
            for trial in  meta[subj]["trials"][group]:
                c3dfile = meta[subj]["trials"][group][trial]["c3dfile"]
                c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                osimkey = c3d_extract(trial, c3dfile, c3dpath, lab, xdir, threshold, ref_model)
    return None



'''
c3d_extract(trial, c3dpath, c3dpath, lab, xdir, threshold, ref_model):
    Extract the motion data from the C3D file to arrays, and returns a dict
    containing all the relevant file metadata, force data and marker data.
'''
def c3d_extract(trial, c3dfile, c3dpath, lab, xdir, threshold, ref_model):
    
    # load C3D file
    itf = c3d.c3dserver()
    c3d.open_c3d(itf, c3dpath + "/" + c3dfile)       
    
    # get all file metadata, and all force plate and video C3D data
    fmeta = c3d.get_dict_groups(itf)
    fforces = c3d.get_dict_forces(itf, frame=True, time=True)
    fmarkers = c3d.get_dict_markers(itf, frame=True, time=True)
    
    # subject name
    fname= fmeta["SUBJECTS"]["NAMES"][0]
    if not fname: fname = "NoName"
    
    # C3D key with all data from C3D file
    c3dkey = C3DKey(fname, fmeta, fforces, fmarkers)

    # trial data only from C3D key
    trialkey = TrialKey(lab, c3dkey, xdir)
    
    # opensim data
    osimkey = OpenSimKey(trialkey, ref_model, c3dpath, threshold)
    
    # save key files
    with open(os.path.join(c3dpath, trial + "_c3dkey.pkl"),"wb") as f: pk.dump(c3dkey, f)
    with open(os.path.join(c3dpath, trial + "_trialkey.pkl"),"wb") as g: pk.dump(trialkey, g)
    with open(os.path.join(c3dpath, trial + "_osimkey.pkl"),"wb") as h: pk.dump(osimkey, h)
    
    return osimkey
    



'''
change_coordinates(oldvec, rotmat, originvec):
    Perform coordinate transformation and change of origin
'''
def change_coordinates(oldvec, rotmat, originvec):
    return originvec + np.matmul(rotmat, oldvec)



'''
find_fp_centre_from_lab_origin(corners):
    Find force plate centre in lab (Vicon) coordinates
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
calculate_centre_of_pressure_fp(n, vc2o, F, M):
    Calculate centre of pressure in force plate coordinates using equations
    from Tim Dorn's GaitExtract toolbox. Note: the GaitExtract documentation
    contains an error in the equations on page 14. The correct equations are
    given in the Matlab function getKinetics.m
        
        Given:
        
            (cop_x, cop_y, cop_z): vector from origin to centre of pressure
            (a, b, c): vector from geometric centre to origin
            (F_x, F_y, F_z): force plate force
            (M_x, M_y, M_z): force plate moment
                
        Equations:
            
            cop_x = -1 * ( (M_y + c * F_x) / F_z ) + a
            cop_y = ( (M_x - c * F_y) / F_z ) + b
            cop_z = 0 (default)        
'''
def calculate_centre_of_pressure_fp(ns, vc2o, F, M):
    cop = np.zeros([ns,3])
    for n in range(ns):
        cop[n,0] = -1 * ((M[n,1] + vc2o[2] * F[n,0]) / F[n,2]) + vc2o[0]
        cop[n,1] = ((M[n,0] - vc2o[2] * F[n,1]) / F[n,2]) + vc2o[1]
        if np.isnan(cop[n,0]) or np.isnan(cop[n,1]):
            cop[n,0] = 0
            cop[n,1] = 0
    return cop



'''
calculate_vertical_free_moment(ns, vc2o, F, M, cop):
    Calculate vertical free moment in force plate coordinates using equations
    from Tim Dorn's GaitExtract toolbox.
        
        Given:
        
            (T_x, T_y, T_z): vector of free moment
            (cop_x, cop_y, cop_z): vector from origin to centre of pressure
            (a, b, c): vector from geometric centre to origin
            (F_x, F_y, F_z): force plate force
            (M_x, M_y, M_z): force plate moment
                
        Equations:
            
            T_x = 0 (default)
            T_y = 0 (default)
            T_z = M_z - ((cop_x - a) * F_y) + ((cop_y - b) * F_x)       
'''
def calculate_vertical_free_moment(ns, vc2o, F, M, cop):
    T = np.zeros([ns,3])
    for n in range(ns):
        T[n,2] = M[n,2] - ((cop[n,0] - vc2o[0]) * F[n,1]) + ((cop[n,1] - vc2o[1]) * F[n,0])
    return T

