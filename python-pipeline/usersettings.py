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
TRAILSettings:
    UserSettings for LASEM TRAIL project.
'''
class TRAILSettings(UserSettings):
    def __init__(self):
        
        self.isempty = False
        
        
        
        # ******************************
        # GENERAL SETTINGS
        
        # data folders
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data"
        self.infolder = r"inputDatabase\Events"
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["Baseline"]

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_opensim_results_all_"
        
        # meta data file
        self.metadatafile = "TRAIL.pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run_stance"] = ["EP", "FAST"]
        self.trialprefixes["run_stridecycle"] = ["EP", "FAST"]
        self.trialprefixes["run_stridecycle_fast"] = ["FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"
        
        # output samples
        self.samples = 101
       
        
       
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
        self.fom_scalefactor = 1.5
        self.lst_scalefactor = -1
        
        # OpenSim RRA parameters
        self.rraiter = 2
                
        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.5   # minimum -0.03 sec
        self.cmc_end_time_offset = -0.05     # due to final event time errors
        