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
        
        # file prefixes
        self.subjprefix = []
        self.trialprefix = []
        
        
        
        
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
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\lasem-trail-running\python-pipeline\reference_model"
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        self.trialprefix = ""
        