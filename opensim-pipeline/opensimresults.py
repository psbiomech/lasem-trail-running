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
    Data storage class containing all OpenSim output data. Data is resampled to
    the desired number of samples.
'''
class OsimResultsKey():
    def __init__(self, osimkey, trialkey, analyses, user, nsamp):
        self.subject = osimkey.subject
        self.trial = osimkey.trial
        self.age = osimkey.age
        self.mass = osimkey.mass
        self.model = osimkey.model
        self.lab = osimkey.lab
        self.task = osimkey.task
        self.dataset = osimkey.dataset
        self.condition = osimkey.condition
        self.events = osimkey.events
        self.outpath = osimkey.outpath
        self.nsamp = nsamp
        self.__calc_discrete_spatiotemporal(osimkey, trialkey, user)
        self.__get_results_raw(osimkey, analyses, nsamp)
        self.__get_results_split(osimkey, analyses, user, nsamp)
        return None        

    def __calc_discrete_spatiotemporal(self, osimkey, trialkey, user):
        
        # Estimate average trial speed (this may only be meaningful for some
        # kinds of tasks, e.g. walking and running)
        avgtrialspeed = np.linalg.norm(osimkey.markers[user.avg_trialspeed_marker][-1, 0:2] - osimkey.markers[user.avg_trialspeed_marker][0, 0:2]) / (osimkey.markers["time"][-1] - osimkey.markers["time"][0])


        # Estimate cadence from step time based on task, only if relevant
        # use events from TrialKey to maximise number of foot-strike events
        
        # Task: RUN
        if self.task.casefold() == "run":
            
            # Get the foot-strikes and indices
            fsidxs = [fn for fn, f in enumerate(trialkey.events["labels"]) if f.endswith("FS")]
            fsfeet = [trialkey.events["labels"][f][0] for f in fsidxs]
            
            # Count foot-strikes not including the first. Ideally this should
            # just be len(fsidx)-1, however, some fast trials are missing 
            # events, so need to count manually.
            nfs = 0
            for f in range(0, len(fsfeet) - 1):
                if ((fsfeet[f] == "R") and (fsfeet[f + 1] == "L")):
                    nfs = nfs + 1
                elif ((fsfeet[f] == "L") and (fsfeet[f + 1] == "R")):
                    nfs = nfs + 1
                elif (fsfeet[f] == fsfeet[f + 1]):      # missing event
                    nfs = nfs + 2
                    
            # Time window between first and last event
            timewindow = trialkey.events["time"][fsidxs[-1]] - trialkey.events["time"][fsidxs[0]]
            avgsteptime = timewindow / nfs
            freq = 1 / avgsteptime
            cadence = freq * 60
            
        # Otherwise, set to nominal value
        else:
            cadence = -1.0
            
                
        self.avgtrialspeed = avgtrialspeed
        self.cadence = cadence
        
        return None
        
    def __get_results_raw(self, osimkey, analyses, nsamp):
        
        # initialise dict
        results = {}
       
        # file suffix and extension
        filext = {}
        filext["ik"] = "_ik.mot"
        filext["id"] = "_id.sto"
        filext["so"] = "_so_force.sto"
        filext["rra"] = []
        filext["cmc"] = []
        filext["jr"] = "_jr_ReactionLoads.sto"
        filext["bk"] = "_bk_pos_global.sto"     # TBD: add vel and acc
        filext["emg"] = "_emg_envelopes.sto"
        filext["grf"] = "_trimmed_grf.sto"
        
        # header rows
        # note: may differ from actual number of header rows as pandas skips
        # some blank initial rows
        headnum = {}
        headnum["ik"] = 8
        headnum["id"] = 6
        headnum["so"] = 10
        headnum["rra"] = []
        headnum["cmc"] = []
        headnum["jr"] = 9
        headnum["bk"] = 13
        headnum["emg"] = 6
        headnum["grf"] = 6
       
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
            
            # Load data into df, if it exists (e.g. not all participants have EMG)
            datafile = os.path.join(osimkey.outpath, ans, osimkey.trial + filext[ans])
            if os.path.isfile(datafile): 
                
                # Raw data
                datadf = pd.read_csv(datafile, sep="\t", header=headnum[ans])
                            
                # Split data and headers
                headers = datadf.columns.tolist()
                data = datadf.to_numpy()
                
                # resample data
                datanew = resample1d(data, nsamp)
            
            # If no data, just return None type
            else:
                datanew = []
                headers = []
            
            # store in dict
            results[ans] = {}
            results[ans]["data"] = datanew
            results[ans]["headers"] = headers
        
        self.results = {}
        self.results["raw"] = results     
            
        return None
    
    def __get_results_split(self, osimkey, analyses, user, nsamp):
        
        # initialise dict
        results = {}
        
        # results processing parameters
        flip = user.results_flip
        columns = user.results_columns
        headers = user.results_headers
        
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
            
            # initialise dict
            results[ans] = {}        
            
            # split by feet
            for f, foot in enumerate(["r", "l"]):                             
                
                results[ans][foot] = {}
                
                # copy raw data
                data0 = None
                data0 = self.results["raw"][ans]["data"].copy()                
                
                # ###################################
                # ADDITIONAL PROCESSING BASED ON TASK AND DATASET
                #
                # Includes left leg trial flipping, as the way this is done can
                # be trial-dependent
                
                # If no data, then return empty list (e.g. not all participants
                # have EMG data)
                if len(data0) == 0:
                    results[ans][foot]["data"] = []
                    results[ans][foot]["headers"] = []
                    results[ans][foot]["events"] = []
                    results[ans][foot]["event_times"] = []
                    results[ans][foot]["event_steps"] = []
                    continue
                
                # Static trial
                elif self.dataset.casefold() == "static":
                    print("Static trial. Nothing to be done.")
                    continue
                    
                # Task: RUN
                elif self.task.casefold() == "run":
                    
                    # foot not used (i.e. no events and no leg_task assigned)
                    if self.events["leg_task"][f] == "not_used":
                        
                        # store empty in dict
                        results[ans][foot] = {}
                        results[ans][foot]["data"] = np.zeros([nsamp, len(headers[ans])])
                        results[ans][foot]["headers"] = headers[ans]
                        results[ans][foot]["event_times"] = -1
                        results[ans][foot]["event_steps"] = -1
                        continue
                    
                    # foot is used (i.e. events exist and a leg_task is assigned)
                    # match dataset and find time window for foot based on leg_task
                    else:
                    
                        # run stridecycle
                        if self.dataset.casefold() == "run_stridecycle":
                            
                            # flip columns for left leg trials
                            if f == 1:
                                data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)
                            
                            # time window depends on leg task and number of events   
                            if self.events["leg_task"][f] == "stridecycle":
                                e0 = self.events["labels"].index(foot.upper() + "FS")
                                if len(self.events["labels"]) == 3:
                                    e1 = e0 + 2
                                else:
                                    e1 = e0 + 4
                                t0 = self.events["time"][e0]
                                t1 = self.events["time"][e1]
                            else:
                                e0 = self.events["labels"].index(foot.upper() + "FS")
                                e1 = self.events["labels"].index(foot.upper() + "FO")
                                t0 = self.events["time"][e0]
                                t1 = self.events["time"][e1]                        
                        
                        
                        # run stance
                        elif self.dataset.casefold() == "run_stance":
                            
                            # flip columns for left leg trials
                            if f == 1:
                                data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)
                                
                            # time window depends on leg task and number of events 
                            if self.events["leg_task"][f] == "stance":
                                e0 = self.events["labels"].index(foot.upper() + "FS")
                                e1 = self.events["labels"].index(foot.upper() + "FO")
                                t0 = self.events["time"][e0]
                                t1 = self.events["time"][e1]     
                                
                        
                # Task: STEP DOWN AND PIVOT
                elif self.task.casefold() == "sdp":
                    
                    # flip columns for left-foot-first left-turning trials
                    if osimkey.events["labels"][0][0].casefold() == "l":
                        data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)                   
                    
                    # time window
                    e0 = 0
                    e1 = 5
                    t0 = self.events["time"][e0]
                    t1 = self.events["time"][e1]    

                    
                # Task: HOP FOR DISTANCE
                elif self.task.casefold() == "hfd":
                    
                    # foot not used (i.e. no events and no leg_task assigned)
                    if self.events["leg_task"][f] == "not_used":
                        
                        # store empty in dict
                        results[ans][foot] = {}
                        results[ans][foot]["data"] = np.zeros([nsamp, len(headers[ans])])
                        results[ans][foot]["headers"] = headers[ans]
                        continue                    
                                        
                    # Foot is used
                    elif self.events["leg_task"][f] == "hfd":
                        
                            # flip columns for left foot trials
                            if osimkey.events["labels"][1][0].casefold() == "l":
                                data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)                   
                            
                            # time window
                            # Note: for hop for distance extract from first generic
                            # event to last generic event.
                            e0 = 0
                            e1 = 3
                            t0 = self.events["time"][e0]
                            t1 = self.events["time"][e1]  
                                
                #
                # ###################################

                # trim columns
                data0 = data0[:, columns[ans][f]].copy()
                
                # trim rows (time window)
                r00 = np.where(data0[:, 0] <= t0)[0]
                if r00.size == 0:
                    r0 = 0
                else:
                    r0 = r00[-1]
                r1 = np.where(data0[:, 0] <= t1)[0][-1]
                data1 = data0[r0:r1 + 1, :]
                
                # resample data, currently uses simple 1D interpolation but
                # need to find a package that emulates Matlab's resample()
                data = resample1d(data1, nsamp)
                            
                # store in dict
                results[ans][foot]["data"] = data
                results[ans][foot]["headers"] = headers[ans]
                results[ans][foot]["events"] = self.events["labels"][e0:e1+1]
                results[ans][foot]["event_times"] = self.events["time"][e0:e1+1]
                results[ans][foot]["event_steps"] = np.round(nsamp * (self.events["time"][e0:e1+1] - self.events["time"][e0]) / (self.events["time"][e1] - self.events["time"][e0]))
        
        self.results["split"] = results        
            
        return None        
        
        

'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
opensim_results_batch_process(meta, analyses, user, nsamp):
    Batch process OpenSim results text files to OsimResultsKeys.
'''
def opensim_results_batch_process(meta, analyses, user, nsamp, restart = -1):
    
    # extract OpenSim data
    osimkey = {}
    failedfiles = []
    startflag = 0
    for subj in meta:
        
        # skip the study info
        if subj.casefold() == "study": continue

    
        # Skip to restart participant, process until last restart participant.
        # Python uses lazy evaluation so combined expressions are efficient.
        if restart != -1:
            if startflag == 1:
                if (type(restart) == tuple) and (subj == restart[1]):
                    startflag = 0            
            elif startflag == 0:
                if (type(restart) == str) and (subj == restart):
                    startflag = 1
                elif (type(restart) == tuple) and (subj == restart[0]):
                    if restart[0] != restart[1]:
                        startflag = 1
                else:
                    continue    
    
        
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                

                #****** TESTING ******
                #if not (trial == "TRAIL296_EP02"): continue;
                #*********************
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                # ignore MVC trials
                ismvc = meta[subj]["trials"][group][trial]["ismvc"]
                if ismvc: continue

                try:
                            
                    # Load the trial OsimKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath, trial + "_osimkey.pkl")
                    with open(pkfile, "rb") as fid:
                        osimkey = pk.load(fid)

                    # Load the TrialKey (for cadence calculation)
                    pkfile = os.path.join(c3dpath, trial + "_trialkey.pkl")
                    with open(pkfile, "rb") as fid:
                        trialkey = pk.load(fid)                    
                        
                    # get the OpenSim results
                    osimresultskey = OsimResultsKey(osimkey, trialkey, analyses, user, nsamp)
                    
                    # save OsimResultsKey to file
                    with open(os.path.join(c3dpath, trial + "_opensim_results.pkl"), "wb") as f:
                        pk.dump(osimresultskey, f)
                                    
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial) 
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)
                          
    print("\n")                
    
    return failedfiles
    
    

'''
export_opensim_results(meta, user, analyses, nsamp):
    Collate OpenSim results into dataframes and export to text for Rstats.
    Note: nsamp should be the same as that used for OsimResultsKey
    otherwise the trial is skipped.
        nsamp: number of samples if different from user.nsamp (int)
        normalise: normalise biomechanical data (bool, default = False)
        emgnormalise: normalise tp "peak", "mvc", or "none" (default)
'''
def export_opensim_results(meta, user, analyses, nsamp, normalise = False, emgnormalise = "none"):
    
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        # skip the study info
        if subj.casefold() == "study": continue        
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        # Rectify subject name format for consistency with TRAIL baseline
        # Current: TRAIL123, Required: TRAIL_123
        subjcorrected = subj[0:5] + "_" + subj[5:8]
        
        # Sex, group and knee reference
        sex = meta[subj]["sex"]
        gtype = meta[subj]["type"]   # Group type: control or surgical (group already used so needed another var name)
        knee = meta[subj]["knee"]
        
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # for MVC normalisation for EMG, load the collated MVC data
            if emgnormalise.casefold() == "mvc":
                mvcpkl = os.path.join(meta[subj]["outpath"], group, subj + "_MVC.pkl")
                with open(mvcpkl,"rb") as fid0:
                    mvcs = pk.load(fid0)            
            
            
            # process dynamic trials only
            for trial in meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                # ignore MVC trials
                ismvc = meta[subj]["trials"][group][trial]["ismvc"]
                if ismvc: continue
            
                # Rectify trial name format for consistency with TRAIL baseline
                # Current: TRAIL123_FAST08, Required: TRAIL_123_FAST08
                trialcorrected = trial[0:5] + "_" + trial[5:]
            
                try:
                
                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                       
                    
                    # trial info
                    task = osimresultskey.task
                    dataset = osimresultskey.dataset
                    condition = osimresultskey.condition
                    avgtrialspeed = osimresultskey.avgtrialspeed
                    cadence = osimresultskey.cadence

                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # leg data window, i.e. leg_task (terrible name)
                        data_type = osimresultskey.events["leg_task"][f]
                        
                        # if the leg has no data, then ignore
                        if data_type.casefold() == "not_used": continue
                        
                        # Check if symptomatic surgical limb or not, or control
                        if gtype=="surgical":
                            if knee=="bilateral":
                                issurgical = 1
                            elif foot==knee[0]:
                                issurgical = 1
                            else:
                                issurgical = 0
                        else:
                            issurgical = -1
                        
                        # analysis
                        for ans in analyses:
                            
                            # ignore scaling
                            if ans.casefold() == "scale": continue
                        
                            # data array
                            data = osimresultskey.results["split"][ans][foot]["data"]
                            varheader = osimresultskey.results["split"][ans][foot]["headers"]
                            
                            # Skip if no data (e.g. not all participants have
                            # EMG data)
                            if len(data) == 0: continue
                            
                            # resample data if necessary
                            if nsamp != osimresultskey.nsamp:
                                data = resample1d(data, nsamp)
                        
                            # Get the events
                            event_times = osimresultskey.results["split"][ans][foot]["event_times"]
                            event_steps = osimresultskey.results["split"][ans][foot]["event_steps"]
                        
                            # Pad events list if required
                            if dataset == "run_stridecycle":
                                if len(event_times) < 5:
                                    event_times = np.concatenate([event_times, [-1]*(5-len(event_times))])
                                    event_steps = np.concatenate([event_steps, [-1]*(5-len(event_steps))])
                        
                        
                            # variable
                            for v, variable in enumerate(varheader):
                                
                                # data for the variable (includes time)
                                drow = data[:, v]

                                # Determine normalisation coefficient
                                # (TBD, set to default 1, as we don't have mass
                                # or height data from TRAIL yet)
                                normfactor = 1
                                if normalise:
                                    if variable.casefold() == "time":
                                        normfactor = 1
                                    elif ans == "ik":
                                        normfactor = 1
                                    elif ans == "id":
                                        if variable.casefold().startswith("pelvis_t"):
                                            normfactor = 1
                                            #normfactor = 1 / mass * user.gravity
                                        else:
                                            normfactor = 1
                                            #normfactor = 100 / (mass * user.gravity * height)
                                    elif ans == "bk":
                                        normfactor = 1
                                    
                                    # Normalise if required
                                    drow = drow  * normfactor
    
    
                                # EMG normalisation
                                if (ans == "emg") and (variable != "time"):
                                    if emgnormalise.casefold() == "peak":
                                        drow = drow / np.max(np.abs(drow))
                                    elif emgnormalise.casefold() == "mvc":
                                        ftvar = foot + variable
                                        mvc = np.abs(mvcs["rms"][ftvar.upper()])
                                        if mvc == -1:    # no MVC data available, just normalise to peak
                                            drow = drow / np.max(np.abs(drow))
                                        else:           # normalise to mean MVC
                                            drow = drow / mvc
                                                
                                # create new line of data
                                csvrow = [subjcorrected, trialcorrected, sex, gtype, task, dataset, condition, data_type, foot, issurgical, avgtrialspeed, cadence, ans, variable] + event_times.tolist() + event_steps.tolist() + drow.tolist()
                                csvdata.append(csvrow)
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "trial", "sex", "group", "task", "dataset", "condition", "data_type", "data_leg", "is_surgical", "avg_speed_m/s", "cadence_steps/min", "analysis", "variable"] + ["etime" + str(n) for n in range(0, len(event_times))] + ["esteps" + str(n) for n in range(0, len(event_steps))] + ["t" + str(n) for n in range(1, nsamp + 1)]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvfileprefix + meta["study"]["task"] + "_" + meta["study"]["dataset"] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, meta["study"]["task"], meta["study"]["dataset"], user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles






'''
export_opensim_results_subject_mean(meta, user, analyses, nsamp):
    Collate OpenSim results into dataframes and export to text for Rstats.
    Note: nsamp should be the same as that used for OsimResultsKey
    otherwise the trial is skipped.
        nsamp: number of samples if different from user.nsamp (int)
        normalise: normalise biomechanical data (bool, default = False)
        emgnormalise: normalise tp "peak", "mvc", or "none" (default)
'''
def export_opensim_results_subject_mean(meta, user, analyses, nsamp, normalise = False, emgnormalise = "none"):
    
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        # skip the study info
        if subj.casefold() == "study": continue        
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        # Rectify subject name format for consistency with TRAIL baseline
        # Current: TRAIL123, Required: TRAIL_123
        subjcorrected = subj[0:5] + "_" + subj[5:8]
        
        # Sex, group and knee reference
        sex = meta[subj]["sex"]
        gtype = meta[subj]["type"]   # Group type: control or surgical (group already used so needed another var name)
        knee = meta[subj]["knee"]
        
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      

            # for MVC normalisation for EMG, load the collated MVC data
            if emgnormalise.casefold() == "mvc":
                mvcpkl = os.path.join(meta[subj]["outpath"], group, subj + "_MVC.pkl")
                with open(mvcpkl,"rb") as fid0:
                    mvcs = pk.load(fid0)  
            
            # process dynamic trials only
            for trial in meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                # ignore MVC trials
                ismvc = meta[subj]["trials"][group][trial]["ismvc"]
                if ismvc: continue
            
                # Rectify trial name format for consistency with TRAIL baseline
                # Current: TRAIL123_FAST08, Required: TRAIL_123_FAST08
                trialcorrected = trial[0:5] + "_" + trial[5:]
            
                try:
                
                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                        
                    # trial info
                    task = osimresultskey.task
                    dataset = osimresultskey.dataset
                    condition = osimresultskey.condition
                    avgtrialspeed = osimresultskey.avgtrialspeed
                    cadence = osimresultskey.cadence

                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # leg data window, i.e. leg_task (terrible name)
                        data_type = osimresultskey.events["leg_task"][f]
                        
                        # if the leg has no data, then ignore
                        if data_type.casefold() == "not_used": continue
                        
                        # Check if symptomatic surgical limb or not, or control
                        if gtype=="surgical":
                            if knee=="bilateral":
                                issurgical = 1
                            elif foot==knee[0]:
                                issurgical = 1
                            else:
                                issurgical = 0
                        else:
                            issurgical = -1
                        
                        # analysis
                        for ans in analyses:
                            
                            # ignore scaling
                            if ans.casefold() == "scale": continue
                        
                            # data array
                            data = osimresultskey.results["split"][ans][foot]["data"]
                            varheader = osimresultskey.results["split"][ans][foot]["headers"]
                            
                            # Skip if no data (e.g. not all participants have
                            # EMG data)
                            if len(data) == 0: continue
                            
                            # resample data if necessary
                            if nsamp != osimresultskey.nsamp:
                                data = resample1d(data, nsamp)
                        
                            # Get the events
                            event_times = osimresultskey.results["split"][ans][foot]["event_times"]
                            event_steps = osimresultskey.results["split"][ans][foot]["event_steps"]
                        
                            # Pad events list if required
                            if dataset == "run_stridecycle":
                                if len(event_times) < 5:
                                    event_times = np.concatenate([event_times, [-1]*(5-len(event_times))])
                                    event_steps = np.concatenate([event_steps, [-1]*(5-len(event_steps))])
                        
                        
                            # variable
                            for v, variable in enumerate(varheader):
                                
                                # data for the variable (includes time)
                                drow = data[:, v]
   
                                # Determine normalisation coefficient
                                # (TBD, set to default 1, as we don't have mass
                                # or height data from TRAIL yet)
                                normfactor = 1
                                if normalise:
                                    if variable.casefold() == "time":
                                        normfactor = 1
                                    elif ans == "ik":
                                        normfactor = 1
                                    elif ans == "id":
                                        if variable.casefold().startswith("pelvis_t"):
                                            normfactor = 1
                                            #normfactor = 1 / mass * user.gravity
                                        else:
                                            normfactor = 1
                                            #normfactor = 100 / (mass * user.gravity * height)
                                    elif ans == "bk":
                                        normfactor = 1
                                    
                                    # Normalise if required
                                    drow = drow  * normfactor
    
    
                                # EMG normalisation
                                if (ans == "emg") and (variable != "time"):
                                    if emgnormalise.casefold() == "peak":
                                        drow = drow / np.max(np.abs(drow))
                                    elif emgnormalise.casefold() == "mvc":
                                        ftvar = foot + variable
                                        mvc = np.abs(mvcs["max"][ftvar.upper()])
                                        if mvc == -1:    # no MVC data available, just normalise to peak
                                            drow = drow / np.max(np.abs(drow))
                                        else:           # normalise to mean MVC
                                            drow = drow / mvc
   
    
                                # create new line of data
                                csvrow = [subjcorrected, trialcorrected, sex, gtype, task, dataset, condition, data_type, foot, issurgical, avgtrialspeed, cadence, ans, variable] + event_times.tolist() + event_steps.tolist() + drow.tolist()
                                csvdata.append(csvrow)
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "trial", "sex", "group", "task", "dataset", "condition", "data_type", "data_leg", "is_surgical", "avg_speed_m/s", "cadence_steps/min", "analysis", "variable"] + ["etime" + str(n) for n in range(0, len(event_times))] + ["esteps" + str(n) for n in range(0, len(event_steps))] + ["t" + str(n) for n in range(1, nsamp + 1)]
    csvdf = pd.DataFrame(csvdata, columns = headers)
    
    # Drop non-numeric columns
    csvdf = csvdf.drop(["trial", "sex", "data_leg"], axis=1)
    
    # Group
    csvdf_grouped = csvdf.groupby(["subject", "group", "task", "dataset", "condition", "data_type", "is_surgical", "analysis", "variable"])

    # Descriptives
    csvdf_grouped_mean = csvdf_grouped.mean().reset_index()
    csvdf_grouped_sd = csvdf_grouped.std().reset_index()

    # Rearrange dataframes to interleave mean and sd rows (much easier in dplyr 
    # with relocate()!)
    csvdf_grouped_mean["statistic"] = "mean"
    dfmean = csvdf_grouped_mean.pop("statistic")
    csvdf_grouped_mean.insert(csvdf_grouped_mean.columns.get_loc("variable") + 1, dfmean.name, dfmean)      
    csvdf_grouped_sd["statistic"] = "sd"
    dfsd = csvdf_grouped_sd.pop("statistic")
    csvdf_grouped_sd.insert(csvdf_grouped_sd.columns.get_loc("variable") + 1, dfsd.name, dfsd)  
    
    # Interleave mean and sd rows
    csvdf_grouped_mean["sortidx"] = range(0, len(csvdf_grouped_mean))
    csvdf_grouped_sd["sortidx"] = range(0, len(csvdf_grouped_sd))
    csvdf_descriptives = pd.concat([csvdf_grouped_mean, csvdf_grouped_sd]).sort_values(["sortidx", "statistic"])
    csvdf_descriptives.drop(columns = "sortidx", inplace = True)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvdescfileprefix + meta["study"]["task"] + "_" + meta["study"]["dataset"] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, meta["study"]["task"], meta["study"]["dataset"], user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf_descriptives.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles





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