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
        
        # data folders
        self.rootpath = []
        self.infolder = []
        self.outfolder = []
        self.trialgroupfolders = []
        self.csvfolder = []
        
        # meta data file
        self.metadatafile = []
        
        # reference model
        self.refmodelpath = []
        self.refmodelfile = []

        # setup file folders
        self.refsetuppath = []
        self.refexternalloads = []
        self.refsetupscale = []
        self.refsetupik = []
        self.refsetupid = []

        # OpenSim analysis codes
        self.scalecode = []
        self.ikcode = []
        self.idcode = []
        
        # file prefixes
        self.subjprefix = []

        # static trial info
        self.staticprefix = []
        self.staticused = []
        self.staticfpchannel = []        
        
        # file suffixes based on task
        self.trialprefixes = {}
        
        # file name format regex pattern
        self.fnpat = []      

        # output samples
        self.samples = []        
    
    
        
'''
TRAILSettings:
    UserSettings for LASEM TRAIL project.
'''
class TRAILSettings(UserSettings):
    def __init__(self):
        
        # data folders
        self.rootpath = r"C:\Users\Owner\Documents\data\TRAIL Test Data"
        self.infolder = r"inputDatabase\Events"
        self.outfolder = "outputDatabase"
        self.trialgroupfolders = ["Baseline"]

        # export data
        self.csvfolder = "csvfolder"
        self.csvfile = "trail_opensim_results_all.csv"
        
        # meta data file
        self.metadatafile = "TRAIL.pkl"
        
        # reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-model"
        self.refmodelfile = "LASEM_TRAIL_ReferenceModel.osim"

        # setup file and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensim-reference-setup"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refsetupscale = "LASEM_TRAIL_Setup_Scale.xml"
        self.refsetupik = "LASEM_TRAIL_Setup_IK.xml"
        self.refsetupid = "LASEM_TRAIL_Setup_ID.xml"
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run"] = ["EP", "FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"
        
        # output samples
        self.samples = 101