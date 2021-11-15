# -*- coding: utf-8 -*-
"""
Lab set ups for C3D extract

@author: Prasanna Sritharan
"""


import numpy as np


'''
---------------------------------
------------ CLASSES ------------
---------------------------------
'''


'''
LabKey:
    Storage class for lab set up info
'''
class LabKey():
    def __init__(self, lab_name, n_fp, n_fp_channels, fp_used, fp_dict_name_prefix, fp_channel_prefixes, fp_channel_suffixes, transform_lab_to_opensim_xdir, transform_fp_to_lab, transform_mat_lab_to_opensim, transform_mat_fp_to_lab, marker_list):
        self.lab_name = lab_name
        self.n_fp = n_fp
        self.n_fp_channels = n_fp_channels
        self.fp_used = fp_used
        self.fp_dict_name_prefix = fp_dict_name_prefix
        self.fp_channel_prefixes = fp_channel_prefixes
        self.fp_channel_suffixes = fp_channel_suffixes
        self.transform_lab_to_opensim_xdir = transform_lab_to_opensim_xdir
        self.transform_fp_to_lab = transform_fp_to_lab
        self.transform_mat_lab_to_opensim = transform_mat_lab_to_opensim
        self.transform_mat_fp_to_lab = transform_mat_fp_to_lab
        self.marker_list = marker_list



'''
-----------------------------------
------ LAB SET-UP FUNCTIONS -------
-----------------------------------
'''


'''
lab_lasem_trail():
    Create a LabKey for the LASEM TRAIL project
'''
def lab_lasem_trail():

    lab_name = "lasem_trail"    

    # force plates
    n_fp = 4;
    n_fp_channels = 6
    fp_used = [3, 4]
    fp_dict_name_prefix = "FP"
    fp_channel_prefixes = ["Force.Fx", "Force.Fy", "Force.Fz", "Moment.Mx", "Moment.My", "Moment.Mz"]
    fp_channel_suffixes = ["", "", "", "", "", ""]
    
    # transform vectors
    transform_lab_to_opensim_xdir = [1, 3, -2]
    transform_fp_to_lab = [-1, 2, -3]
    
    # transform matrices from transform vectors
    transform_mat_lab_to_opensim = create_transform_set_lab_to_opensim(transform_lab_to_opensim_xdir)
    transform_mat_fp_to_lab = create_transform_matrix(transform_fp_to_lab)
    
    # markers
    marker_list = []
    
    # create a lab
    return LabKey(lab_name, n_fp, n_fp_channels, fp_used, fp_dict_name_prefix, fp_channel_prefixes, fp_channel_suffixes, transform_lab_to_opensim_xdir, transform_fp_to_lab, transform_mat_lab_to_opensim, transform_mat_fp_to_lab, marker_list)
    


'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
create_transform_set_lab_to_opensim():
    Create lab to opensim rotation matrices from transform vectors
'''
def create_transform_set_lab_to_opensim(tvec_x):
    
    # initiliase set of transform for directions: X -X Y -Y Z -Z
    # note: Z transforms are assumed and may need to be respecified depending
    #   on the lab set-up
    transform_set = [None]*6
    
    # Vicon X ---> OpenSim X
    transform_set[0] = create_transform_matrix(tvec_x)
    
    # Vicon -X ---> OpenSim X
    tvec_nx = np.multiply(tvec_x,[-1, 1, -1])
    transform_set[1] = create_transform_matrix(tvec_nx)
    
    # Vicon Y ---> Opensim X
    tvec_y = [-tvec_x[2], tvec_x[1], tvec_x[0]]
    transform_set[2] = create_transform_matrix(tvec_y)
    
    # Vicon -Y ---> OpenSim X
    tvec_ny = np.multiply(tvec_y,[-1, 1, -1])
    transform_set[3] = create_transform_matrix(tvec_ny)    
    
    # Vicon Z ---> OpenSim X (assumed, change if necessary)
    tvec_z = [tvec_x[1], -tvec_x[2], -tvec_x[0]]
    transform_set[4] = create_transform_matrix(tvec_z)  

    # Vicon -Z ---> OpenSim X (assumed, change if necessary)
    tvec_nz = np.multiply(tvec_z,[-1, 1, -1])
    transform_set[5] = create_transform_matrix(tvec_nz)
    
    return transform_set


    
'''
create_transform_matrix():
    Create transform matrix from transform vectors
'''
def create_transform_matrix(tvec):
    tmat = np.zeros((3,3))
    for r in range(3):
        tmat[r,abs(tvec[r])-1] = 1*np.sign(tvec[r])   
    return tmat
