# -*- coding: utf-8 -*-
"""
Load and format OpenSim results

@author: Prasanna Sritharan
"""

import os
import pandas as pd
import numpy as np
import pickle as pk


'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''


'''
OsimResultsKey:
    Data storage class containing all OpenSim output data, raw and normalised
    to both BW and %BW*HT.
'''
class OsimResultsKey():
    def __init__(self, osimkey, analyses):
        self.subject = osimkey.subject
        self.trial = osimkey.trial
        self.age = osimkey.age
        self.mass = osimkey.mass
        self.model = osimkey.model
        self.lab = osimkey.lab
        self.task = osimkey.task
        self.outpath = osimkey.outpath
        self.__get_results(osimkey, analyses)
        return None
        
    def __get_results(self, osimkey, analyses):
        
        # initialise dict
        results = {}
       
        # file suffix and extension
        filext = {}
        filext["ik"] = "_ik.mot"
        filext["id"] = "_id.sto"
        filext["so"] = []
        filext["rra"] = []
        filext["cmc"] = []
        filext["jr"] = []
        
        # header rows
        # note: may differ from actual number of header rows as pandas skips
        # some initial rows (there's probably a way to override this, TBD)
        headnum = {}
        headnum["ik"] = 8
        headnum["id"] = 6
        headnum["so"] = []
        headnum["rra"] = []
        headnum["cmc"] = []
        headnum["jr"] = []
       
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
            
            # load data 
            datafile = os.path.join(osimkey.outpath, ans, osimkey.trial + filext[ans])
            datadf = pd.read_csv(datafile, sep="\t", header=headnum[ans])
            headers = datadf.columns.tolist()
            data = datadf.to_numpy()
            
            # store in dict
            results[ans] = {}
            results[ans]["data"] = data
            results[ans]["headers"] = headers
        
        self.results = results
        
            
        return None



'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
opensim_results_batch_process(meta):
    Batch process OpenSim results text files to OsimResultsKeys.
'''
def opensim_results_batch_process(meta, analyses):
    
   # extract OpenSim data
   osimkey = {}
   for subj in meta:

       print("%s" % "*" * 30)
       print("SUBJECT: %s" % subj)
       print("%s" % "*" * 30)
       print("\n")
       
       for group in meta[subj]["trials"]:
           
           print("Group: %s" % group)
           print("%s" % "=" * 30)                      
           
           # process dynamic trials only
           for trial in  meta[subj]["trials"][group]:                
               
               # ignore static trials
               isstatic = meta[subj]["trials"][group][trial]["isstatic"]
               if isstatic: continue
           
               print("Dynamic trial: %s" % trial)
           
               # load the trial OsimKey
               c3dpath = meta[subj]["trials"][group][trial]["outpath"]
               pkfile = os.path.join(c3dpath,trial + "_osimkey.pkl")
               with open(pkfile,"rb") as fid:
                   osimkey = pk.load(fid)
                   
               # get the OpenSim results
               osimresultskey = OsimResultsKey(osimkey, analyses)
               
               # save OsimResultsKey to file
               with open(os.path.join(c3dpath, trial + "_opensim_results.pkl"),"wb") as f:
                   pk.dump(osimresultskey, f)
                   
   
   print("\n")                
   
   return None
    
    

'''
collate_opensim_results(meta):
    Collate OpenSim results into dataframes and export to text for Rstats.
'''
def collate_opensim_results(meta, analyses):
    pass