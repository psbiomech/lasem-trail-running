# -*- coding: utf-8 -*-
"""
Plot participant waveforms for data checking - knee angle and moment

@author: Prasanna Sritharan
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pk
import usersettings as uset

# user settings
user = uset.TRAILSettings_RUN()

# data file
srcpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\run\run_stance\csvfolder"
srcfile = "trail_opensim_results_ikid_run_run_stance.csv"

# output file
outpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\\run\run_stance\datacheck"
if not os.path.isdir(outpath): os.makedirs(outpath)



# %% PREPARE DATA

# load meta data
dbfilepath = os.path.join(user.rootpath, user.outfolder, "run", "run_stance", user.metadatafile)
with open(dbfilepath, "rb") as fid:
    meta = pk.load(fid)

# load CSV data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# retain only knee angle and knee moment rows
df = df0[(df0["variable"] == "knee_angle") | (df0["variable"] == "knee_angle_moment")]



# %% PLOT OUTPUT


print("Generating plots for data check...")

# task: run
task = "run"

# dataset: run_stridecycle, run_stance
dataset = "run_stance"

# data leg: r (ignore left leg)
leg = "r"

# condition: ep, fast
for c in ["ep", "fast"]:
    
    # data_type: stridecycle, stance
    for d in ["stridecycle", "stance"]:
              
        # variable: knee_angle, knee_angle_moment
        for v in ["knee_angle", "knee_angle_moment"]:
            
            print("\n----------------------------------------")
            print("PARAMETERS: %s, %s, %s, %s" % (task.upper(), c.upper(), d.upper(), v.upper()))
            print("----------------------------------------")            
            
            # subjects
            for subj in meta:
                
                # skip the study info
                if subj.casefold() == "study": continue                
                
                print("%s" % subj, end = "")
                
                # output folder
                figpath = os.path.join(outpath, task, subj)
                if not os.path.isdir(figpath): os.makedirs(figpath)
                
                # get the participant data for the given parameter combination
                data = df[(df["data_leg"] == leg) & (df["condition"] == c) & (df["data_type"] == d) & (df["variable"] == v) & (df["subject"] == subj)]
                
                # drop info columns except trial name, transpose
                dnames = data.loc[:, "trial"]
                data = data.drop(columns = ["subject", "trial", "task", "dataset", "condition", "data_type", "data_leg", "analysis", "variable"])
                data = data.T
                data.columns = dnames
                      
                # if no data, skip
                if data.empty:
                    print(" *** NO DATA ***")
                    continue
                else:
                    print("")
                
                # plot participant data
                data.plot(title = subj + "_" + v + "_" + d + "_" + c + "_" + leg)
                
                # save figure
                figfile = subj + "_" + v + "_" + d + "_" + c + "_" + leg + ".pdf"
                plt.savefig(os.path.join(figpath, figfile))
                
                # close figure
                plt.close();


print("\nPlot generation complete.\n")