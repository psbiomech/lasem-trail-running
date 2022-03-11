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
        #self.metadatafile = self.project + ".pkl"     
        
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
        self.filter_cutoff = 40
        self.filter_threshold = 16   # prefer 15N, but 16N required for EP_08
        self.smooth_cop_fixed_offset = 50
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
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints.xml"
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 2.0
        self.fom_scalefactor["iliacus"] = 3.0
        self.fom_scalefactor["psoas"] = 3.0
        self.fom_scalefactor["ehl"] = 3.0
        self.fom_scalefactor["edl"] = 3.0
        self.fom_scalefactor["vas"] = 3.0
        self.fom_scalefactor["glm"] = 3.0
        self.fom_scalefactor["bf"] = 3.0
        self.fom_scalefactor["semi"] = 3.0       
        self.fom_scalefactor["gas"] = 1.0
        self.fom_scalefactor["sol"] = 5.0
        self.lst_scalefactor = {}
        self.lst_scalefactor["all"] = 0      
        
        # OpenSim RRA parameters
        self.rraiter = 2
        self.rra_end_time_offset = 0.03
                
        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation



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
        #self.metadatafile = self.project + ".pkl"        
        
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
        self.filter_cutoff = 40
        self.filter_threshold = 16
        self.smooth_cop_fixed_offset = 50
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
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints.xml"
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 2.0
        self.fom_scalefactor["iliacus"] = 3.0
        self.fom_scalefactor["psoas"] = 3.0
        self.fom_scalefactor["ehl"] = 3.0
        self.fom_scalefactor["edl"] = 3.0
        self.fom_scalefactor["vas"] = 3.0
        self.fom_scalefactor["glm"] = 3.0
        self.fom_scalefactor["bf"] = 3.0
        self.fom_scalefactor["semi"] = 3.0       
        self.fom_scalefactor["gas"] = 1.5
        self.fom_scalefactor["sol"] = 6.0
        self.lst_scalefactor = {}
        self.lst_scalefactor["all"] = 0      
        
        # OpenSim RRA parameters
        self.rraiter = 2       
        self.rra_end_time_offset = 0.03        
        
        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation