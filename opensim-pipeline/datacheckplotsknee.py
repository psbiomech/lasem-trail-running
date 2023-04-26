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
user = uset.TRAILSettings_RUN_ID()

# data file
srcpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\csvfolder"
srcfile = "trail_run_opensim_results_ikid.csv"

# output file
outpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase\datacheck"
if not os.path.isdir(outpath): os.makedirs(outpath)



# %% PREPARE DATA

# load meta data
dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    meta = pk.load(fid)

# load CSV data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# retain only knee angle and knee moment rows
df = df0[(df0["variable"] == "knee_angle") | (df0["variable"] == "knee_angle_moment")]



# %% PLOT OUTPUT


print("Generating plots for data check...")

# task: run_stridecycle, run_stance
task = "run_stance"

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
                
                print("%s" % subj, end = "")
                
                # output folder
                figpath = os.path.join(outpath, task, subj)
                if not os.path.isdir(figpath): os.makedirs(figpath)
                
                # get the participant data for the given parameter combination
                data = df[(df["data_leg"] == leg) & (df["condition"] == c) & (df["data_type"] == d) & (df["variable"] == v) & (df["subject"] == subj)]
                
                # drop info columns except trial name, transpose
                dnames = data.loc[:, "trial"]
                data = data.drop(columns = ["subject", "trial", "task", "condition", "data_type", "data_leg", "analysis", "variable"])
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


# # empty output list of lists
# # (create the output table as a list of lists, then convert to dataframe
# # as iteratively appending new dataframe rows is computationally expensive)
# csvdata = []
    
# # extract OpenSim data
# print("Collating data into lists...\n")
# failedfiles = []
# for subj in meta:

#     print("%s" % "*" * 30)
#     print("SUBJECT: %s" % subj)
#     print("%s" % "*" * 30)
            
#     for group in meta[subj]["trials"]:
        
#         print("Group: %s" % group)
#         print("%s" % "=" * 30)                      
        
#         # process dynamic trials only
#         for trial in  meta[subj]["trials"][group]:                
            
#             # ignore static trials
#             isstatic = meta[subj]["trials"][group][trial]["isstatic"]
#             if isstatic: continue
        
#             try:
            
#                 # load the trial OsimResultsKey
#                 c3dpath = meta[subj]["trials"][group][trial]["outpath"]
#                 pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
#                 with open(pkfile,"rb") as fid:
#                     osimresultskey = pk.load(fid)
                    
#                 # trial task
#                 task = user.results_task_for_output
                
#                 # condition
#                 condition = osimresultskey.condition

#                 # foot
#                 for f, foot in enumerate(["r","l"]):
                    
#                     # leg data window, i.e. leg_task (terrible name)
#                     data_type = osimresultskey.events["leg_task"][f]
                    
#                     # analysis
#                     for ans in analyses:
                        
#                         # ignore scaling
#                         if ans.casefold() == "scale": continue
                    
#                         # data array
#                         data = osimresultskey.results["split"][ans][foot]["data"]
#                         varheader = osimresultskey.results["split"][ans][foot]["headers"]
                    
#                         # variable
#                         for v, variable in enumerate(varheader):
                            
#                             # data for the variable (includes time)
#                             drow = data[:, v]

#                             # create new line of data
#                             csvrow = [subj, trial, task, condition, data_type, foot, ans, variable] + drow.tolist()
#                             csvdata.append(csvrow)
            
#             except:
#                 print("Dynamic trial: %s *** FAILED ***" % trial)
#                 failedfiles.append(trial)
#             else:
#                 print("Dynamic trial: %s" % trial)

# # create dataframe
# print("\nCreating dataframe...")
# headers = ["subject", "trial", "task", "condition", "data_type", "data_leg", "analysis", "variable"] + ["t" + str(n) for n in range(1,102)]
# csvdf = pd.DataFrame(csvdata, columns = headers)

# # write data to file with headers
# print("\nWriting to CSV text file...")
# csvfile = user.csvfileprefix + ".csv"
# fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
# if not os.path.exists(fpath): os.makedirs(fpath)
# csvdf.to_csv(os.path.join(fpath,csvfile), index = False)


# print("\n")