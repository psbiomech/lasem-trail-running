# -*- coding: utf-8 -*-
"""
User settings parameters

@author: Prasanna Sritharan
"""



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''


'''
UserSettings:
    Template class for user settings for C3D processing pipeline.
'''
class UserSettings():
    def __init__(self):
        
        self.isempty = True    
    
    
        
'''
TRAILSettings_RUN_EP(UserSettings):
    UserSettings for LASEM TRAIL Project: RUN EP
'''
class TRAILSettings_RUN_EP(UserSettings):
    def __init__(self):
        
        self.isempty = False
        
                
        # ******************************
        # GENERAL SETTINGS

        # project
        self.project = "TRAIL_RUN_EP"
        
        # data folders
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data"
        self.infolder = r"inputDatabase\Events"
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["Baseline"]

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_run_ep_opensim_results_all_"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"     
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run_stance"] = ["EP"]
        self.trialprefixes["run_stridecycle"] = ["EP"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     
       
        # force plate data filter
        self.filter_butter_order = 4
        self.filter_cutoff = 15
        self.filter_threshold = 16   # prefer 15N, but 16N required for EP_08
        self.smooth_cop_fixed_offset = 50    # required but not used
        self.smooth_window = 25

       
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline"
        self.logfile = "out.log"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-model"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel.osim"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-setup"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
        # OpenSim additional files
        self.additionalfilesfolder = "EP"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN_EP.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN_EP.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN_EP.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN_EP.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN_EP.xml"
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 2.5
        self.fom_scalefactor["sol"] = 4.0
        self.lom_scalefactor = {}
        self.lom_scalefactor["all"] = 1.1  
        
        # OpenSim IK parameters
        self.kinematics_filter_cutoff = 15
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2
        self.rra_start_time_offset = -0.05
        self.rra_end_time_offset = 0.05
        self.prescribe_upper_body_motion = True
        self.prescribed_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.0



'''
TRAILSettings_RUN_FAST(UserSettings):
    UserSettings for LASEM TRAIL Project: RUN FAST
'''
class TRAILSettings_RUN_FAST(UserSettings):
    def __init__(self):
        
        self.isempty = False
        
                
        # ******************************
        # GENERAL SETTINGS
        
        # project
        self.project = "TRAIL_RUN_FAST"
        
        # data folders
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data"
        self.infolder = r"inputDatabase\Events"
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["Baseline"]

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_run_fast_opensim_results_all_"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run_stance"] = ["FAST"]
        self.trialprefixes["run_stridecycle"] = ["FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     
       
        # force plate data filter
        self.filter_butter_order = 4
        self.filter_cutoff = 15
        self.filter_threshold = 15
        self.smooth_cop_fixed_offset = 25   # required but not currently used
        self.smooth_window = 20
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline"
        self.logfile = "out.log"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-model"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel_Unclamped_NoUpperActuators.osim"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-setup"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
        # OpenSim additional files
        self.additionalfilesfolder = "FAST"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN_FAST.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN_FAST.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN_FAST.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN_FAST.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN_FAST.xml"
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 3.0
        #self.fom_scalefactor["sol"] = 5.0        
        #self.fom_scalefactor["vas"] = 5.0
        #self.lom_scalefactor = {}
        #self.lom_scalefactor["sol"] = 1.2

        
        # OpenSim IK parameters
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribed_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.0