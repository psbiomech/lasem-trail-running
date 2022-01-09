# -*- coding: utf-8 -*-
"""
Load and format OpenSim results

@author: Prasanna Sritharan
"""

import os
import pandas as pd
import numpy as np
import pickle as pk
from scipy.interpolate import interp1d


'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''


'''
OsimResultsKey:
    Data storage class containing all OpenSim output data, raw and normalised
    to both BW and %BW*HT. Data is resample to the desired number of samples.
'''
class OsimResultsKey():
    def __init__(self, osimkey, analyses, nsamp):
        self.subject = osimkey.subject
        self.trial = osimkey.trial
        self.age = osimkey.age
        self.mass = osimkey.mass
        self.model = osimkey.model
        self.lab = osimkey.lab
        self.task = osimkey.task
        self.outpath = osimkey.outpath
        self.__get_results(osimkey, analyses, nsamp)
        return None
        
    def __get_results(self, osimkey, analyses, nsamp):
        
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
            
            # resample data
            datanew = resample1d(data, nsamp)
            
            # store in dict
            results[ans] = {}
            results[ans]["data"] = datanew
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



'''
resample1d(data, nsamp):
    Simple resampling by 1-D interpolation (rows = samples, cols = variable).
    Data can be a 1-D or multiple variables in a 2D array-like object.
'''
def resample1d(data, nsamp):

    # convert list to array
    if isinstance(data, list):
        data = np.array([data]).transpose()
        ny = 1        

    # data dimensions
    nx = data.shape[0]
    ny = data.shape[1]

    # old sample points
    x = np.linspace(0, nx - 1, nx)
        
    # new sample points
    xnew = np.linspace(0, nx - 1, nsamp)
    
    # parse columns
    datanew = np.zeros([nsamp, ny])
    for col in range(0, ny):
        
        # old data points
        y = data[:, col]
        
        # convert to cubic spline function
        fy = interp1d(x, y, kind = "cubic", fill_value = "extrapolate")
    
        # new data points
        ynew = fy(xnew)
    
        # store column
        datanew[:, col] = ynew
        
    
    return datanew