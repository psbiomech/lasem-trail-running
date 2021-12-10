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
        
        # model folders
        self.refmodelfile = []
        self.refmodelpath = []

        # setup file folders
        self.refsetuppath = []
        self.refsetupscale = []
        
        # file prefixes
        self.subjprefix = []

        # static file
        self.staticfile = []        
        
        # file suffixes based on task
        self.trialprefixes = {}
        
        # file name format regex pattern
        self.fnpat = []        
        
    
    
        
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
        
        # model folders
        self.refmodelfile = "Rajagopal2015.osim"
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensom-reference-model"

        # setup file folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\opensom-reference-setup"
        self.refsetupscale = "ACLR_Reference_Setup_Scale.xml"
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static file
        self.staticfile = "Static_01"        
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run"] = ["EP", "FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"