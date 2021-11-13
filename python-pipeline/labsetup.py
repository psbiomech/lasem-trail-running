# -*- coding: utf-8 -*-
"""
Lab set ups for C3D extract

@author: Prasanna Sritharan
"""


'''
LabKey:
    Storage class for lab set up info
'''
class LabKey():
    def __init__(self, nfp, fpused, fpprefixes, fpsuffixes, lab_to_opensim, fp_to_lab, marker_list):
        self.nfp = nfp
        self.fpused = fpused
        self.fpprefixes = fpprefixes
        self.fpsuffixes = fpsuffixes
        self.transform_lab_to_opensim = lab_to_opensim
        self.transform_fp_to_lab = fp_to_lab
        self.marker_list = marker_list




'''
LabKey: LASEM TRAIL Project
    Storage class for lab set up info
'''
def lab_lasem_trail():

    # force plates
    nfp = 4;
    used_force_plates = [3,4]
    force_prefixes = ["Force.Fx","Force.Fy","Force.Fz","Moment.Mx","Moment.My","Moment.Mz"],
    force_suffixes = ["","","","","",""]
    
    # coordinate systems
    lab_to_opensim = [1, 3, -2]
    fp_to_lab = [-1, 2, -3]
    
    # markers
    marker_list = []
    
    # create a lab
    return LabKey(nfp, used_force_plates, force_prefixes, force_suffixes, lab_to_opensim, fp_to_lab, marker_list)
    