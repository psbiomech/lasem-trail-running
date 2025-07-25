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
    Parent user settings class for C3D processing pipeline.
'''
class UserSettings():
    def __init__(self):
        
        
        # ******************************
        # GENERAL SETTINGS        
        
        # Add generic settings here
  


        # ******************************
        # C3D DATA PROCESSING  
        
        # Add generic parameters here
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        self.jrcode = "jr"
        self.bkcode = "bk"
        self.emgcode = "emg"
        
        # limb code
        self.leg = ["r", "l"]
        
        
        
'''
TRAILSettings_RUN(UserSettings):
    UserSettings for LASEM TRAIL Project: RUN
'''
class TRAILSettings_RUN(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(TRAILSettings_RUN, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL"
        self.infolder = ["inputDatabase"]
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["BASELINE"]  
        
        # project
        self.project = "TRAIL_RUN"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_opensim_results_"
        self.csvdescfileprefix = "trail_opensim_results_subject_descriptives_"
        
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # MVC prefix
        self.mvcprefix = "MVC"
        
        # C3D file suffixes for datasets based on task: RUN
        self.trialprefixes = {}
        self.trialprefixes["run"] = {}
        self.trialprefixes["run"]["run_stance"] = ["EP", "FAST"]
        self.trialprefixes["run"]["run_stridecycle"] = ["EP", "FAST"]
        self.trialprefixes["run"]["run_stance_ep"] = ["EP"]
        self.trialprefixes["run"]["run_stridecycle_ep"] = ["EP"]
        self.trialprefixes["run"]["run_stance_fast"] = ["FAST"]
        self.trialprefixes["run"]["run_stridecycle_fast"] = ["FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL\d+_(EP|FAST|STATIC)\d+"
        self.tasktoknum = 1   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     

        # Analog channels
        self.analogchannelnames = {}
                
        # EMG subcohort list file
        self.emglistfile = "EMG tracking.xlsx"
              
        # marker data filter (set cutoff to -1 if not required)
        self.marker_filter_butter_order = 4
        self.marker_filter_cutoff = -1
       
        # force plate data filter (set cutoff to -1 if not required)
        self.fp_filter_butter_order = 4
        self.fp_filter_cutoff = 15
        self.fp_smooth_transitions = False
        self.fp_filter_threshold = -1
        self.fp_smooth_cop_fixed_offset = 0   # required but not currently used
        self.fp_expand_window = 20
                        
        # marker to use for estimating trial speed
        self.avg_trialspeed_marker = "SACR"
        
        
        # ******************************
        # OPENSIM PARAMETERS

        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline"
        self.triallogfolder = "log"
        self.logfile = "opensim.log"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-setup"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-models"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel_Unclamped.osim"
        self.refmodeltposefile = "LASEM_TRAIL_ReferenceModel_Tpose_Unclamped.osim"
        self.tposestaticlistfile = "TRAIL T-pose static trials.xlsx"
        
        # OpenSim additional files
        self.additionalfilesfolder = "RUN"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN.xml"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 3.0
        #self.lom_lmt_scalefactor = {}
        #self.lom_lmt_scalefactor["all"] = 1.1
        
        # OpenSim IK parameters
        self.kinematics_filter_butter_order = 4
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.03
        
        # OpenSim JR parameters
        self.jr_joints = {}
        self.jr_joints["all"] = ["child", "child"]
        self.jr_use_cmc_forces = False   
        
        # OpenSim BK parameters
        self.bk_bodies = ["all"]
        self.bk_output_com = True
        self.bk_output_in_local_frame = False
        self.bk_use_cmc_results = False
        
        
        
        # ******************************
        # OPENSIM RESULTS
        
        # task name for output
        self.results_task_for_output = "run"
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 28, 29]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = []  
        self.results_flip["bk"] = [3, 4, 5, 45, 46, 47, 51, 52, 53, 57, 58, 59, 63, 64, 65, 69, 70, 71, 75, 76, 77, 81, 82, 83, 111, 112, 113, 117, 118, 119, 123, 124, 125, 129, 130, 131, 135]
        self.results_flip["emg"] = []
        self.results_flip["grf"] = [9, 12, 16, 17]        
        
        # foot columns (incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = []   
        self.results_columns["bk"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 79, 80, 81, 82, 83, 84],
                                      [0, 1, 2, 3, 4, 5, 6, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]]
        self.results_columns["emg"] = [[0, 1, 2, 3, 4, 5, 6, 7],
                                       [0, 8, 9, 10, 11, 12, 13, 14]]
        self.results_columns["grf"] = [[0, 1, 2, 3, 4, 5, 6, 13, 14, 15],
                                       [0, 7, 8, 9, 10, 11, 12, 16, 17, 18]]
                
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = []  
        self.results_headers["bk"] = ["time", "pelvis_X", "pelvis_Y", "pelvis_Z", "pelvis_Ox", "pelvis_Oy", "pelvis_Oz", "femur_X", "femur_Y", "femur_Z", "femur_Ox", "femur_Oy", "femur_Oz", "tibia_X", "tibia_Y", "tibia_Z", "tibia_Ox", "tibia_Oy", "tibia_Oz", "patella_X", "patella_Y", "patella_Z", "patella_Ox", "patella_Oy", "patella_Oz", "talus_X", "talus_Y", "talus_Z", "talus_Ox", "talus_Oy", "talus_Oz", "calcn_X", "calcn_Y", "calcn_Z", "calcn_Ox", "calcn_Oy", "calcn_Oz", "toes_X", "toes_Y", "toes_Z", "toes_Ox", "toes_Oy", "toes_Oz", "torso_X", "torso_Y", "torso_Z", "torso_Ox", "torso_Oy", "torso_Oz"]
        self.results_headers["emg"] = ["time", "sol", "gaslat", "gasmed", "semiten", "bflh", "vasmed", "vaslat"]
        self.results_headers["grf"] = ["time", "grf_vx", "grf_vy", "grf_vz", "cop_px", "cop_py", "cop_pz", "grm_mx", "grf_my", "grm_mz"]



        # ******************************
        # ADDITIONAL ANALYSES
        
        # outfile file prefixes
        self.csvfileprefix_analyses_jap = "trail_analyses_results_jap"
        self.csvfileprefix_analyses_jaw = "trail_analyses_results_jaw"





'''
TRAILSettings_RUN_EMG(UserSettings):
    UserSettings for LASEM TRAIL Project: RUN (EMG subcohort)
'''
class TRAILSettings_RUN_EMG(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(TRAILSettings_RUN_EMG, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL"
        self.infolder = ["inputDatabase"]
        self.outfolder = "outputDatabase_emg"
        self.trialgroupfolders = ["BASELINE"]  
        
        # project
        self.project = "TRAIL_RUN"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_opensim_results_"
        self.csvdescfileprefix = "trail_opensim_results_subject_descriptives_"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # MVC prefix
        self.mvcprefix = "MVC"
        
        # C3D file suffixes for datasets based on task: RUN
        self.trialprefixes = {}
        self.trialprefixes["run"] = {}
        self.trialprefixes["run"]["run_stance"] = ["EP", "FAST"]
        self.trialprefixes["run"]["run_stridecycle"] = ["EP", "FAST"]
        self.trialprefixes["run"]["run_stance_ep"] = ["EP"]
        self.trialprefixes["run"]["run_stridecycle_ep"] = ["EP"]
        self.trialprefixes["run"]["run_stance_fast"] = ["FAST"]
        self.trialprefixes["run"]["run_stridecycle_fast"] = ["FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL\d+_(EP|FAST|STATIC|MVC)\w+"
        self.tasktoknum = 1   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101
       
        
         
        # ******************************
        # C3D DATA PROCESSING     
        
        # Analog channels
        self.analogchannelnames = {"RSOL": "Sensor 10.EMG10",
                                   "RGASLAT": "Sensor 4.EMG4",
                                   "RGASMED": "Sensor 9.EMG9",
                                   "RSEMITEN": "Sensor 5.EMG5",
                                   "RBFLH": "Sensor 6.EMG6",
                                   "RVASMED": "Sensor 13.EMG13",
                                   "RVASLAT": "Sensor 14.EMG14",                                   
                                   "LSOL": "Sensor 8.EMG8",
                                   "LGASLAT": "Sensor 16.EMG16",
                                   "LGASMED": "Sensor 7.EMG7",
                                   "LSEMITEN": "Sensor 1.EMG1",
                                   "LBFLH": "Sensor 2.EMG2",
                                   "LVASMED": "Sensor 11.EMG11",
                                   "LVASLAT": "Sensor 12.EMG12"}
                
        # EMG subcohort list file
        self.emglistfile = "EMG tracking.xlsx"
        
        # EMG data filter  (set cutoff to -1 if not required)
        # Set to -1 if running EMG processing in OpenSim pipeline
        self.emg_filter_butter_order = 4
        self.emg_filter_cutoff = 10
        
        # Normalise EMG data ("none", "peak", "mvc")
        self.emg_normalise = "none"
        
        # EMG use Hilbert transform for envelope
        # Set to False if running EMG processing in OpenSim pipeline
        self.emg_use_hilbert = False
 
        
 
        # MVC groupings
        self.mvcgroupings = {"RIGHTCALF": ["RSOL", "RGASMED", "RGASLAT"],
                             "RIGHTHAMS": ["RSEMITEN", "RBFLH"],
                             "RIGHTQUAD": ["RVASMED", "RVASLAT"],
                             "LEFTCALF": ["LSOL", "LGASMED", "LGASLAT"],
                             "LEFTHAMS": ["LSEMITEN", "LBFLH"],
                             "LEFTQUAD": ["LVASMED", "LVASLAT"]}
    
        self.mvcsamplewindow = [0.35, 0.65]   # proportion of trial
        self.mvcnsamp = 10001
         
         
        
        # marker data filter (set cutoff to -1 if not required)
        self.marker_filter_butter_order = 4
        self.marker_filter_cutoff = -1
       
        # force plate data filter (set cutoff to -1 if not required)
        self.fp_filter_butter_order = 4
        self.fp_filter_cutoff = 15
        self.fp_smooth_transitions = False
        self.fp_filter_threshold = -1
        self.fp_smooth_cop_fixed_offset = 0   # required but not currently used
        self.fp_expand_window = 20
                        
        # marker to use for estimating trial speed
        self.avg_trialspeed_marker = "SACR"
        
        
        # ******************************
        # OPENSIM PARAMETERS

        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline"
        self.triallogfolder = "log"
        self.logfile = "opensim.log"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-setup"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-models"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel_Unclamped.osim"
        self.refmodeltposefile = "LASEM_TRAIL_ReferenceModel_Tpose_Unclamped.osim"
        self.tposestaticlistfile = "TRAIL T-pose static trials.xlsx"
        
        # OpenSim additional files
        self.additionalfilesfolder = "RUN"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN.xml"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 3.0
        #self.lom_lmt_scalefactor = {}
        #self.lom_lmt_scalefactor["all"] = 1.1
        
        # OpenSim IK parameters
        self.kinematics_filter_butter_order = 4
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.03
        
        # OpenSim JR parameters
        self.jr_joints = {}
        self.jr_joints["all"] = ["child", "child"]
        self.jr_use_cmc_forces = False   
        
        # OpenSim BK parameters
        self.bk_bodies = ["all"]
        self.bk_output_com = True
        self.bk_output_in_local_frame = False
        self.bk_use_cmc_results = False
        
        
        # EMG processing parameters
        self.emg_process_envelope = "movingrms"   # moving average
        self.emg_process_convolve_window = 200   # samples
        
        
        
        
        
        
        # ******************************
        # OPENSIM RESULTS
        
        # task name for output
        self.results_task_for_output = "run"
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 28, 29]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = []  
        self.results_flip["bk"] = [3, 4, 5, 45, 46, 47, 51, 52, 53, 57, 58, 59, 63, 64, 65, 69, 70, 71, 75, 76, 77, 81, 82, 83, 111, 112, 113, 117, 118, 119, 123, 124, 125, 129, 130, 131, 135]
        self.results_flip["emg"] = []
        self.results_flip["grf"] = [9, 12, 16, 17]        
        
        # foot columns (incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = []   
        self.results_columns["bk"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 79, 80, 81, 82, 83, 84],
                                      [0, 1, 2, 3, 4, 5, 6, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]]
        self.results_columns["emg"] = [[0, 1, 2, 3, 4, 5, 6, 7],
                                       [0, 8, 9, 10, 11, 12, 13, 14]]
        self.results_columns["grf"] = [[0, 1, 2, 3, 4, 5, 6, 13, 14, 15],
                                       [0, 7, 8, 9, 10, 11, 12, 16, 17, 18]]
                
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = []  
        self.results_headers["bk"] = ["time", "pelvis_X", "pelvis_Y", "pelvis_Z", "pelvis_Ox", "pelvis_Oy", "pelvis_Oz", "femur_X", "femur_Y", "femur_Z", "femur_Ox", "femur_Oy", "femur_Oz", "tibia_X", "tibia_Y", "tibia_Z", "tibia_Ox", "tibia_Oy", "tibia_Oz", "patella_X", "patella_Y", "patella_Z", "patella_Ox", "patella_Oy", "patella_Oz", "talus_X", "talus_Y", "talus_Z", "talus_Ox", "talus_Oy", "talus_Oz", "calcn_X", "calcn_Y", "calcn_Z", "calcn_Ox", "calcn_Oy", "calcn_Oz", "toes_X", "toes_Y", "toes_Z", "toes_Ox", "toes_Oy", "toes_Oz", "torso_X", "torso_Y", "torso_Z", "torso_Ox", "torso_Oy", "torso_Oz"]
        self.results_headers["emg"] = ["time", "sol", "gaslat", "gasmed", "semiten", "bflh", "vasmed", "vaslat"]
        self.results_headers["grf"] = ["time", "grf_vx", "grf_vy", "grf_vz", "cop_px", "cop_py", "cop_pz", "grm_mx", "grf_my", "grm_mz"]




        # ******************************
        # ADDITIONAL ANALYSES
        
        # outfile file prefixes
        self.csvfileprefix_analyses_jap = "trail_analyses_results_jap"
        self.csvfileprefix_analyses_jaw = "trail_analyses_results_jaw"









        
        
'''
TRAILSettings_HFD(UserSettings):
    UserSettings for LASEM TRAIL Project: HOP FOR DISTANCE
'''
class TRAILSettings_HFD(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(TRAILSettings_HFD, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL"
        self.infolder = ["inputDatabase"]
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["BASELINE"]  
        
        # project
        self.project = "TRAIL_HFD"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_opensim_hfd_results_ikid"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # C3D file suffixes for datasets based on task: RUN
        self.trialprefixes = {}
        self.trialprefixes["hfd"] = {}
        self.trialprefixes["hfd"]["hfd"] = ["HFD_RIGHT", "HFD_LEFT"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL\d+_(HFD_RIGHT|HFD_LEFT|STATIC)\d+"
        self.tasktoknum = 1   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     

        # Analog channels
        self.analogchannelnames = {}
                
        # EMG subcohort list file
        self.emglistfile = "EMG tracking.xlsx"
              
        # marker data filter (set cutoff to -1 if not required)
        self.marker_filter_butter_order = 4
        self.marker_filter_cutoff = -1
       
        # force plate data filter (set cutoff to -1 if not required)
        self.fp_filter_butter_order = 4
        self.fp_filter_cutoff = 15
        self.fp_smooth_transitions = False
        self.fp_filter_threshold = -1
        self.fp_smooth_cop_fixed_offset = 0   # required but not currently used
        self.fp_expand_window = 20
        
        # marker to use for estimating trial speed
        self.avg_trialspeed_marker = "SACR"
        
        
        # ******************************
        # OPENSIM PARAMETERS

        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline"
        self.triallogfolder = "log"
        self.logfile = "opensim.log"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-setup"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-models"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel_Unclamped.osim"
        
        # OpenSim additional files
        self.additionalfilesfolder = "HFD"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_HFD.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_HFD.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_HFD.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_HFD.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_HFD.xml"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 3.0
        #self.lom_lmt_scalefactor = {}
        #self.lom_lmt_scalefactor["all"] = 1.1
        
        # OpenSim IK parameters
        self.kinematics_filter_butter_order = 4
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.03
        
        # OpenSim JR parameters
        self.jr_joints = {}
        self.jr_joints["all"] = ["child", "child"]
        self.jr_use_cmc_forces = False        
        
        
        # ******************************
        # OPENSIM RESULTS
        
        # task name for output
        self.results_task_for_output = "hfd"
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 28, 29]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = []  
                
        # foot columns (incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = []   
        
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = []           
        
        
        
# '''
#  TRAILSettings_RUN_ID(UserSettings):
#      UserSettings for LASEM TRAIL Project: RUN using ID model
# '''
# class TRAILSettings_RUN_ID(UserSettings):
#     def __init__(self):
        
        
#         # inherit parent attributes
#         super(TRAILSettings_RUN_ID, self).__init__()
        
                        
#         # ******************************
#         # GENERAL SETTINGS
        
#         self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL"
#         self.infolder = ["inputDatabase"]
#         self.outfolder = "outputDatabase"
#         self.trialgroupfolders = ["BASELINE"]  
        
#         # project
#         self.project = "TRAIL_RUN"

#         # export data
#         self.csvfolder = "csvfolder"
#         self.csvfileprefix = "trail_run_opensim_results_ikid"
        
#         # meta data file
#         self.metadatafile = self.project + ".pkl"        
        
#         # file prefixes
#         self.subjprefix = "TRAIL"
        
#         # static trial info
#         self.staticprefix = "STATIC"
#         self.staticused = "Static_01"
#         self.staticfpchannel = "Force.Fz3"
        
#         # file suffixes based on task
#         self.trialprefixes = {}
#         self.trialprefixes["run_stance"] = ["EP", "FAST"]
#         self.trialprefixes["run_stridecycle"] = ["EP", "FAST"]
#         self.trialprefixes["run_stance_ep"] = ["EP"]
#         self.trialprefixes["run_stridecycle_ep"] = ["EP"]
#         self.trialprefixes["run_stance_fast"] = ["FAST"]
#         self.trialprefixes["run_stridecycle_fast"] = ["FAST"]
                        
#         # file name format regex pattern:
#         #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
#         self.fnpat = "TRAIL\d+_(EP|FAST|STATIC)\d+"
#         self.tasktoknum = 1   # the token + 1 that represents the task name/type
        
#         # output samples
#         self.samples = 101
       
        
#         # ******************************
#         # C3D DATA PROCESSING     
       
#         # marker data filter (set cutoff to -1 if not required)
#         self.marker_filter_butter_order = 4
#         self.marker_filter_cutoff = -1
       
#         # force plate data filter (set cutoff to -1 if not required)
#         self.fp_filter_butter_order = 4
#         self.fp_filter_cutoff = 15
#         self.fp_filter_threshold = -1
#         self.fp_smooth_transitions = False
#         self.fp_smooth_cop_fixed_offset = 0   # required but not currently used
#         self.fp_expand_window = 20
        
        
#         # ******************************
#         # OPENSIM PARAMETERS

#         # OpenSim log file
#         self.logfilepath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline"
#         self.triallogfolder = "log"
#         self.logfile = "opensim.log"

#         # OpenSim setup files and folders
#         self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-setup"
#         self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
#         self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
#         self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
#         self.refsetupso = "LASEM_TRAIL_Setup_Analysis.xml"
#         self.refsetuprra = "LASEM_TRAIL_Setup_RRA.xml"
#         self.refsetupcmc = "LASEM_TRAIL_Setup_CMC.xml"
        
#         # OpenSim reference model
#         #   LASEM_TRAIL_ReferenceModel_InverseDynamics_Final: 
#         #       knee/ankle rotations and adductions clamped, others unclamped
#         self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\opensim-pipeline\opensim-models"
#         self.refmodelfile = "LASEM_TRAIL_ReferenceModel_InverseDynamics_Final.osim"
         
#         # OpenSim additional files
#         self.additionalfilesfolder = "RUN"
#         self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
#         self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
#         self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN.xml"
#         self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN.xml"
#         self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN.xml"
#         self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN.xml"
#         self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN.xml"
        
#         # OpenSim Scale parameters
#         #self.fom_scalefactor = {}
#         #self.fom_scalefactor["all"] = 3.0
#         #self.lom_lmt_scalefactor = {}
#         #self.lom_lmt_scalefactor["all"] = 1.1
        
#         # OpenSim IK parameters
#         self.kinematics_filter_butter_order = 4
#         self.kinematics_filter_cutoff = 6
        
#         # OpenSim RRA parameters
#         self.update_mass = True
#         self.rraiter = 2   
#         self.rra_start_time_offset = -0.03  # to enable CMC initalisation
#         self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
#         self.prescribe_upper_body_motion = True
#         self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

#         # OpenSim CMC parameters
#         self.use_rra_model = True
#         self.use_rra_kinematics = True
#         self.use_fast_target = True
#         self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
#         self.cmc_end_time_offset = 0.03
        
#         # OpenSim JR parameters
#         self.jr_joints = {}
#         self.jr_joints["all"] = ["child", "child"]
#         self.jr_use_cmc_forces = False        
        
        
#         # ******************************
#         # OPENSIM RESULTS
        
#         # task name for output
#         self.results_task_for_output = "run"
        
#         # left leg flip columns (incl. time)
#         self.results_flip = {}
#         self.results_flip["ik"] = [2, 3, 6, 28, 29]
#         self.results_flip["id"] = [2, 3, 6, 14, 15]
#         self.results_flip["so"] = []
#         self.results_flip["rra"] = []
#         self.results_flip["cmc"] = []
#         self.results_flip["jr"] = []  
                
#         # foot columns (incl. time): R, L
#         self.results_columns = {}
#         self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36], 
#                                       [0, 1, 2, 3, 4, 5, 6, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 37, 38, 39, 40, 41, 42, 43]]
#         self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 24, 25, 26, 30, 31, 32, 36, 38, 40, 41], 
#                                       [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 20, 21, 22, 27, 28, 29, 33, 34, 35, 37, 39, 42, 43]]
#         self.results_columns["so"] = []
#         self.results_columns["rra"] = []
#         self.results_columns["cmc"] = []
#         self.results_columns["jr"] = []   
        
#         # headers
#         self.results_headers = {}
#         self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_adduction", "knee_rotation", "ankle_angle", "ankle_adduction", "ankle_rotation", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
#         self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_adduction_moment", "knee_rotation_moment", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "ankle_adduction_moment", "ankle_rotation_moment", "elbow_flex_moment", "pro_sup_moment", "wrist_flex_moment", "wrist_dev_moment"]
#         self.results_headers["so"] = []
#         self.results_headers["rra"] = []
#         self.results_headers["cmc"] = []
#         self.results_headers["jr"] = []