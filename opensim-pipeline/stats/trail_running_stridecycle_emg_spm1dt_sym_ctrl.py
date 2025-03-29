# -*- coding: utf-8 -*-
"""
LASEM RUNNING SPM{t}: STRIDECYCLE

@author: Prasanna Sritharan, August 2023
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
import pickle as pk

# Matplotlib: write text as text not path
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font", **{'family':'sans-serif','sans-serif':['Arial']})


# Data file
srcpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase_emg\run\run_stridecycle\csvfolder"
srcfile = "trail_opensim_results_subject_descriptives_run_run_stridecycle.csv"

# Output file
outpath = r"C:\Users\Owner\Documents\data\TRAIL\outputDatabase_emg\run\run_stridecycle\stats"
if not os.path.isdir(outpath): os.makedirs(outpath)
outfilename = "trail_run_run_stridecycle_emg"

# Exclusions (must have reason)
exclusions = []


# %% PREPARE DATA

# Load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# Remove stdev rows
df1 = df0[df0["statistic"] == "mean"]

# Data labels
analyses = ["emg"]
osimvars = {}
osimvars["emg"] = ["sol", "gaslat", "gasmed", "semiten", "bflh", "vasmed", "vaslat"]
#osimvars["ik"] = ["hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "ankle_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation"]
#osimvars["id"] = [k + "_moment" for k in osimvars["ik"]]
issurgical = [1, -1]
subjtype = ["sym", "ctrl"]
data_type = ["stridecycle"] #["stridecycle", "stance"]
condition = ["ep", "fast"]
#condition =[["stance_more", "stance"], ["stance_more", "stance"]]   # symptomatic stance leg and asymptompatic swing leg

# Get data into arrays
datamat = {}
#eventmat = []
for n in condition:
    datamat[n] = {}
    for d in data_type:
        datamat[n][d] = {}
        df = df1[(df1["condition"]==n) & (df1["data_type"]==d)]
        for a in analyses:
            datamat[n][d][a] = {}
            for v in osimvars[a]:
                
                datamat[n][d][a][v] = {}
                
                # Group 1: sym
                gp1data = df[~df["subject"].isin(exclusions) & (df["is_surgical"] == issurgical[0]) & (df["analysis"] == a) & (df["variable"] == v)]
                datamat[n][d][a][v][0] = gp1data.loc[:, "t1":"t101"].to_numpy()
                       
                # Group 2: ctrl
                gp2data = df[~df["subject"].isin(exclusions) & (df["is_surgical"] == issurgical[1]) & (df["analysis"] == a) & (df["variable"] == v)]
                datamat[n][d][a][v][1] = gp2data.loc[:, "t1":"t101"].to_numpy()
    


# %% EVENTS

# Event time steps and descriptives
# events = {}
# events["data"] = {}
# events["desc"] = {}
# descmat = np.zeros(2)
# for g in range(len(subjtype)):
#     events["data"][g] = eventmat[g]
#     events["desc"][g] = {}
#     events["desc"][g]["mean"] = np.round(np.mean(events["data"][g], axis=0))
#     events["desc"][g]["sd"] = np.round(np.std(events["data"][g], axis=0))
#     descmat[g] = events["desc"][g]["mean"]
# events["desc"]["total"] = {}
# events["desc"]["total"]["mean"] = np.mean(descmat, axis=0)
# events["desc"]["total"]["sd"] = np.std(descmat, axis=0)



# %% RUN ANALYSES: DESCRIPTIVES, SPM{t}

# Calculate group ensemble descriptives from file
desc = {}
for n in condition:
    desc[n] = {}
    for d in data_type:
        desc[n][d] = {}
        for a in analyses:        
            desc[n][d][a] = {}
            for v in osimvars[a]:
                desc[n][d][a][v] = {}
                for s in range(len(subjtype)):    
                    desc[n][d][a][v][s] = {}                      
                    desc[n][d][a][v][s]["mean"] = np.mean(datamat[n][d][a][v][s], axis = 0)
                    desc[n][d][a][v][s]["sd"] = np.std(datamat[n][d][a][v][s], axis = 0)


# F-test: Run SPM{t} and inference for three samples (one-way ANOVA)
bonferroni = 0  # 0=no, 1=yes
significance = [0.05, 0.0167]  # [treat variables as independent vs Bonferroni corrected for 3 comparisons]
spmt = {}
spmtinf = {}
for n in condition:
    spmt[n] = {}
    spmtinf[n] = {}   
    for d in data_type:
        spmt[n][d] = {}
        spmtinf[n][d] = {} 
        for a in analyses: 
            spmt[n][d][a] = {}
            spmtinf[n][d][a] = {}
            for v in osimvars[a]:
                #print("%s, %s, %s, %s" % (n, d, a, v))
                Y0 = datamat[n][d][a][v][0]
                Y1 = datamat[n][d][a][v][1]
                spmt[n][d][a][v] = spm1d.stats.ttest2(Y0, Y1, equal_var=True)
                spmtinf[n][d][a][v] = spmt[n][d][a][v].inference(alpha = significance[bonferroni], two_tailed=True, interp=True)


# Combine for output
sldj = {}
sldj["desc"] = desc
#sldj["events"] = events
sldj["spmt"] = spmt
sldj["spmtinf"] = spmtinf
            
# Pickle it
#with open(os.path.join(outpath, outfilename + ".pkl"),"wb") as f: pk.dump(sldj, f)
    

# %% PLOT OUTPUT

# Plot parameters
# eventlabels = ["PFO1", "PFS2", "NFO1", "NFS2", "PFO3", "PFS4"]
# eventlabelalign = ["left", "right", "left", "right", "left", "right"]
# eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]   

# Figure headers
plotheads = {}
#plotheads["ik"] = ["Hip flexion", "Hip adduction", "Hip rotation", "Knee flexion", "Ankle dorsiflexion", "Lumbar extension", "Lumbar bending", "Lumbar rotation"]
#plotheads["id"] = [k + " moment" for k in plotheads["ik"]]
plotheads["emg"] = ["sol", "gaslat", "gasmed", "semiten", "bflh", "vasmed", "vaslat"]
plotsubjtype = ["symptomatic", "control"]
plotdatatype = ["stride cycle", "stance"]



# Generate plots
#eventlist = events["desc"]["total"]["mean"]
for n in condition:

    for nd, d in enumerate(data_type):    
        
        for a in osimvars:
            
            # Scenario parameters
            #nsubjs = [np.size(datamat[n]["ik"]["time"][s], axis=0) for s in range(len(subjtype))]
    
            # Create plot area
            fig = plt.figure(constrained_layout=True, figsize=(28, 5))   
            fig.suptitle("TRAIL %s: %s (%s)" % (a.upper(), n.upper(), plotdatatype[nd].upper()), fontsize=20)
            spec = fig.add_gridspec(nrows = 2, ncols = len(osimvars["emg"]), height_ratios = [2, 1]) 
            
            # Create plots
            x = range(101)
            for col, v in enumerate(osimvars[a]):                
                
                # Mean
                m0 = desc[n][d][a][v][0]["mean"]
                m1 = desc[n][d][a][v][1]["mean"]
                
                # Upper
                u0 = m0 + desc[n][d][a][v][0]["sd"]
                u1 = m1 + desc[n][d][a][v][1]["sd"]
                
                # Lower
                l0 = m0 - desc[n][d][a][v][0]["sd"]
                l1 = m1 - desc[n][d][a][v][1]["sd"]  
                
                # Plot
                ax = fig.add_subplot(spec[0, col])
                ax.set_title(plotheads[a][col].upper(), fontsize = 12)
                if col == 0:
                    ax.set_ylabel("Normalised activation", fontsize = 12) 
                ax.fill_between(x, l1, u1, alpha = 0.3, linewidth = 0.0, color = "green")
                ax.fill_between(x, l0, u0, alpha = 0.3, linewidth = 0.0, color = "red") 
                ax.plot(x, m1, label = plotsubjtype[1], linewidth = 2.0, color = "green") 
                ax.plot(x, m0, label = plotsubjtype[0], linewidth = 2.0, color = "red")
                ax.set_xlim([x[0], x[-1]])
                ax.set_xlabel("%% of %s" % (plotdatatype[nd]), fontsize = 12)
                if col == 0: ax.legend(frameon = False, loc = "upper right")
                
                # SPM plot
                ax = fig.add_subplot(spec[1, col])
                ax.set_xlabel("%% of %s" % (plotdatatype[nd]), fontsize = 12)
                if col == 0: ax.set_ylabel("SPM{t}", fontsize = 12) 
                spmtinf[n][d][a][v].plot(plot_ylabel = False)
                    
                    
    # save to pdf
    plt.savefig(os.path.join(outpath, outfilename + "_" + n + "_spm1dt_sym_ctrl.pdf"))