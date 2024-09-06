# -*- coding: utf-8 -*-
"""
OpenSim post-hoc analyses: Angular work and power

@author: Prasanna Sritharan
"""



import os
import numpy as np
import pickle as pk
import pandas as pd
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import interp1d
from scipy import signal



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES




'''
-----------------------------------
------- FUNCTIONS: ANALYSES -------
-----------------------------------
'''


'''
analyses_batch_process(meta, user, analyses, filterqs, filter_order, filter_cutoff, restart):
    Batch process post-hoc analyses.
'''
def analyses_batch_process(meta, user, analyses, filterqs = False, filter_order = -1, filter_cutoff = -1, restart = -1):
    
    # extract OpenSim data
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


        print("\n")
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                #****** TESTING ******
                #if not (trial == "TRAIL449_FAST05"): continue;
                #*********************                  
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                try:

                    # Load OsimResultsKey
                    pklpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        osimresultskey = pk.load(fid0)
                            
                    # Output analysis dict
                    analysesdict = {}
                    analysesdict["subject"] = subj
                    analysesdict["group"] = group
                    analysesdict["trial"] = trial
                    analysesdict["outpath"] = pklpath
                                  
                    print("Dynamic trial: %s" % trial)

                    
                    for ans in analyses:
                        
                        # Joint angular impulse
                        if ans == "jai":
                            print("---> Joint angular impulse")
                            analysesdict["joint_angular_impulse"] = calculate_joint_angular_impulse(osimresultskey, user)     
                            
                        # Joint angular power
                        elif ans == "jap":
                            print("---> Joint angular power")
                            analysesdict["joint_angular_power"] = calculate_joint_angular_power(osimresultskey, user, filterqs, filter_order, filter_cutoff)
                            
                        # Joint angular work
                        elif ans == "jaw":
                            print("---> Joint angular work")
                            analysesdict["joint_angular_work"] = calculate_joint_angular_work(osimresultskey, user)                            
                            

                    # Pickle analysis dict
                    with open(os.path.join(pklpath, trial + "_analyses_results.pkl"), "wb") as fid1:
                        pk.dump(analysesdict, fid1)   

                except:
                    print("*** FAILED ***")
                    failedfiles.append(trial)   
                    #raise
                        
                          
    print("\n")                
    
    return failedfiles    



'''
calculate_joint_angular_impulse(osimresultskey, user):
    Calculate the joint angular impulse for each time window specified by
    events, report positive and negative separately.
'''
def calculate_joint_angular_impulse(osimresultskey, user):
    
    # get event times
    events = osimresultskey.events["time"]
    
    # get joint moments
    Tdata = osimresultskey.results["split"][user.idcode]
    
    # calculate joint angular impulse on each leg
    impl = {}
    for leg in user.leg:
        
        # data
        time = Tdata[leg]["data"][:, 0]
        T = Tdata[leg]["data"][:, 1:]
        headers = Tdata[leg]["headers"][1:]
        
        # net
        impl[leg] = {}
        impl[leg]["net"] = {}
        impl[leg]["net"]["headers"] = headers
        impl[leg]["net"]["net"] = np.trapz(T, time, axis = 0)
        
        # net positive
        Tpos = T.copy()
        Tpos[Tpos < 0] = 0
        impl[leg]["net"]["pos"] = np.trapz(Tpos, time, axis = 0)
        
        # net negative
        Tneg = T.copy()
        Tneg[Tneg > 0] = 0
        impl[leg]["net"]["neg"] = np.trapz(Tneg, time, axis = 0)  
        
        # events
        impl[leg]["windows"] = {}
        impl[leg]["windows"]["headers"] = headers
        for e in range(0, len(events) - 1):
            
            # indices for time window
            idx0 = np.where(time >= events[e])[0][0]
            idx1 = np.where(time <= events[e + 1])[0][-1]
            
            # event window data
            timewin = time[idx0:idx1]
            Twin = T[idx0:idx1, :]
            
            # net window
            wlabel = "w" + str(e + 1)
            impl[leg]["windows"][wlabel] = {}
            impl[leg]["windows"][wlabel]["net"] = np.trapz(Twin, timewin, axis = 0)
            
            # net window positive
            Twinpos = Tpos[idx0:idx1, :]
            impl[leg]["windows"][wlabel]["pos"] = np.trapz(Twinpos, timewin, axis = 0)
            
            # net window negative
            Twinneg = Tneg[idx0:idx1, :]
            impl[leg]["windows"][wlabel]["neg"] = np.trapz(Twinneg, timewin, axis = 0)
    
    return impl
    

    
'''
calculate_joint_angular_power(osimresultskey, user, filterqs, filter_order, filter_cutoff):
    Calculate the joint angular power for all available joints using joint angles
    from inverse kinematics and joint moments from inverse dynamics. Filter the
    kinematics if required.
'''
def calculate_joint_angular_power(osimresultskey, user, filterqs = False, filter_order = -1, filter_cutoff = -1):
    
    # Get ID data
    Tdata = osimresultskey.results["split"][user.idcode].copy()
    
    # Get IK data
    Qdata = osimresultskey.results["split"][user.ikcode].copy()
    
    # Calculate joint angular power on each leg
    jap = {}
    for f, leg in enumerate(user.leg):
        
        # Leg joint moments (Nm), remove pelvis forces
        timeT = Tdata[leg]["data"][:, 0]
        T = Tdata[leg]["data"][:, 1:]
        headersT = Tdata[leg]["headers"][1:]
        delidxT = headersT.index("pelvis_tx_force")
        del headersT[delidxT:delidxT+3]
        T = np.delete(T, range(delidxT, delidxT + 3), axis = 1)       
 
        # Leg joint angles (rad), remove pelvis translations
        timeQ = Qdata[leg]["data"][:, 0]
        Q = np.deg2rad(Qdata[leg]["data"][:, 1:])
        headersQ = Qdata[leg]["headers"][1:]
        delidxQ = headersQ.index("pelvis_tx")
        del headersQ[delidxQ:delidxQ+3]
        Q = np.delete(Q, range(delidxQ, delidxQ + 3), axis = 1)   
 
    
        # Initialise output power matrix
        P = np.zeros(np.shape(T))
              
        # Calculate power if the leg is used
        if not osimresultskey.events["leg_task"][f].casefold() == "not_used":
            
            # Filter the joint angles
            Qfilt = Q
            if filterqs:
                sample_rate = 1 / (timeQ[1] - timeQ[0])
                Qfilt = filter_timeseries(Q, sample_rate, filter_order, filter_cutoff)
            
            # Calculate joint velocities (rad/s)
            Qdot0 = np.gradient(Qfilt, timeQ, axis = 0)

            # Filter the joint velocities
            Qdotfilt0 = Qdot0
            if filterqs:
                sample_rate = 1 / (timeQ[1] - timeQ[0])
                Qdotfilt0 = filter_timeseries(Qdot0, sample_rate, filter_order, filter_cutoff)
            
            # Interpolation for joint velocities to match time stamps with moments.
            # Typically the IK time window is expanded to cater for RRA but the ID
            # window is not, so the time stamps do not automatically match.
            Qdot_interpfun = interpolate.interp1d(timeQ, Qdotfilt0, kind = "cubic", axis = 0, fill_value = "extrapolate")
            Qdot = Qdot_interpfun(timeT)
            
            # Calculate joint angular power for each coordinate
            for idxQ, qcoord in enumerate(headersQ):
                
                # Find moment column index of current coordinate as moments and
                # velocities may be out of order
                idxT = [1 if qcoord in head else 0 for i, head in enumerate(headersT)].index(1)
                
                # Get the data for the current coordinate
                qdot = Qdot[:, idxQ]
                t = T[:, idxT]
    
                # Calculate joint angular power: P = Tw (Nm/s)
                P[:, idxQ] = np.multiply(t, qdot) 


        # Add time column to data matrix
        P = np.concatenate((np.reshape(timeQ, (-1, 1)), P), axis = 1)
        
        # Results dict
        jap[leg] = {}
        jap[leg]["data"] = P
        jap[leg]["headers"] = ["time"] + headersQ
    
    return jap



'''
calculate_joint_angular_work(osimresultskey, user):
    Calculate work from power. Output positive and negative work. Waveform
    numerical integration undertaken using Simpson's rule.
'''
def calculate_joint_angular_work(osimresultskey, user):
    
    # Calculate joint angular power
    jap = calculate_joint_angular_power(osimresultskey, user)
    
    # Calculate the joint angular work by integrating waveforms wrt time
    jaw = {}
    for leg in user.leg:
        
        # Split time column and data columns
        t = jap[leg]["data"][:, 0]
        P = jap[leg]["data"][:, 1:]
        
        # NOTE: Numerical integration undertaken using Simpson's rule
        
        # Net work (Nm)
        W = integrate.simpson(P, t, axis = 0)
        
        # Positive work (Nm)
        Ppos = P.copy()
        Ppos[Ppos < 0.0] = 0.0
        Wpos = integrate.simpson(Ppos, t, axis = 0)
        
        # Negative work (Nm)
        Pneg = P.copy()
        Pneg[Pneg > 0.0] = 0.0
        Wneg = integrate.simpson(Pneg, t, axis = 0)
               
        # Results dict
        jaw[leg] = {}
        jaw[leg]["net"] = W
        jaw[leg]["pos"] = Wpos
        jaw[leg]["neg"] = Wneg
        jaw[leg]["headers"] = jap[leg]["headers"][1:]
    
    return jaw





  
  
'''
-----------------------------------
-------- FUNCTIONS: EXPORT --------
-----------------------------------
'''




'''
export_joint_angular_work(meta, user):
    Write joint angular power to Excel.
'''
def export_joint_angular_work(meta, user):
        
    # Empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # Extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        # Skip the study info
        if subj.casefold() == "study": continue        
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        # Rectify subject name format for consistency with TRAIL baseline
        # Current: TRAIL123, Required: TRAIL_123
        subjcorrected = subj[0:5] + "_" + subj[5:8]
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # Process dynamic trials only
            for trial in meta[subj]["trials"][group]:                
 
                #****** TESTING ******
                #if not (trial == "TRAIL001_EP02"): continue;
                #*********************     
                
                # Ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                # Rectify trial name format for consistency with TRAIL baseline
                # Current: TRAIL123_FAST08, Required: TRAIL_123_FAST08
                trialcorrected = trial[0:5] + "_" + trial[5:]
            
                try:
                
                    # Load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
 
                    # Load the trial analyses dict
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_analyses_results.pkl")
                    with open(pkfile,"rb") as fid:
                        analysesdict = pk.load(fid)                       
 
                    # Trial info
                    task = osimresultskey.task
                    dataset = osimresultskey.dataset
                    condition = osimresultskey.condition

                    # Foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # Leg data window, i.e. leg_task (terrible name)
                        data_type = osimresultskey.events["leg_task"][f]
                        
                        # If the leg has no data, then ignore
                        if data_type.casefold() == "not_used": continue
                    
                        # Data array
                        net = analysesdict["joint_angular_work"][foot]["net"]
                        pos = analysesdict["joint_angular_work"][foot]["pos"]
                        neg = analysesdict["joint_angular_work"][foot]["neg"]
                        varheader = analysesdict["joint_angular_work"][foot]["headers"]
                    
                        # Variable
                        for v, variable in enumerate(varheader):
                            
                            # Data for the variable (includes time)
                            drow = [net[v], pos[v], neg[v]]

                            # Create new line of data
                            csvrow = [subjcorrected, trialcorrected, task, dataset, condition, data_type, foot, variable] + drow
                            csvdata.append(csvrow)
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "trial", "task", "dataset", "condition", "data_type", "data_leg", "variable", "net", "positive", "negative"]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvfileprefix_analyses_jaw + "_" + meta["study"]["task"] + "_" + meta["study"]["dataset"] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, meta["study"]["task"], meta["study"]["dataset"], user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles  




'''
export_joint_angular_power(meta, user, nsamp):
    Write joint angular power to Excel.
'''
def export_joint_angular_power(meta, user, nsamp):
        
    # Empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # Extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        # Skip the study info
        if subj.casefold() == "study": continue        
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        # Rectify subject name format for consistency with TRAIL baseline
        # Current: TRAIL123, Required: TRAIL_123
        subjcorrected = subj[0:5] + "_" + subj[5:8]
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # Process dynamic trials only
            for trial in meta[subj]["trials"][group]:                
 
                #****** TESTING ******
                #if not (trial == "TRAIL001_EP02"): continue;
                #*********************     
                
                # Ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                # Rectify trial name format for consistency with TRAIL baseline
                # Current: TRAIL123_FAST08, Required: TRAIL_123_FAST08
                trialcorrected = trial[0:5] + "_" + trial[5:]
            
                try:
                
                    # Load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
 
                    # Load the trial analyses dict
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_analyses_results.pkl")
                    with open(pkfile,"rb") as fid:
                        analysesdict = pk.load(fid)                       
 
                    # Trial info
                    task = osimresultskey.task
                    dataset = osimresultskey.dataset
                    condition = osimresultskey.condition

                    # Foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # Leg data window, i.e. leg_task (terrible name)
                        data_type = osimresultskey.events["leg_task"][f]
                        
                        # If the leg has no data, then ignore
                        if data_type.casefold() == "not_used": continue
                    
                        # Data array
                        data = analysesdict["joint_angular_power"][foot]["data"]
                        varheader = analysesdict["joint_angular_power"][foot]["headers"]
                    
                        # Resample data if necessary
                        # Temporary: disable for run data as results key does not
                        # yet contain nsamp field. Needs to be re-run.
                        # if nsamp != osimresultskey.nsamp:
                        #     data = resample1d(data, nsamp)
                    
                        # Variable
                        for v, variable in enumerate(varheader):
                            
                            # Data for the variable (includes time)
                            drow = data[:, v]

                            # Create new line of data
                            csvrow = [subjcorrected, trialcorrected, task, dataset, condition, data_type, foot, variable] + drow.tolist()
                            csvdata.append(csvrow)
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "trial", "task", "dataset", "condition", "data_type", "data_leg", "variable"] + ["t" + str(n) for n in range(1, nsamp + 1)]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvfileprefix_analyses_jap + "_" + meta["study"]["task"] + "_" + meta["study"]["dataset"] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, meta["study"]["task"], meta["study"]["dataset"], user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles  




'''
export_joint_angular_impulse(meta, user):
    Write joint angular impulses to Excel.
'''
def export_joint_angular_impulse(meta, user):
        
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []   
    
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        print("\n")
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                try:

                    # load results key
                    pklpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        osimresultskey = pk.load(fid0)
                        
                        
                    # movement task and condition
                    movement = osimresultskey.task
                    
                    # windows
                    windows = ["w" + str(e) for e in range(1, len(osimresultskey.events["labels"]))]
                    
                    # pivot leg
                    pivot_leg = osimresultskey.events["labels"][0][0].lower()
                    
                     
                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # pivot leg or non-pivot leg data
                        if foot == pivot_leg:
                            data_leg_role = "pivot"
                        else:
                            data_leg_role = "nonpivot"
                        
                        # joint angular impulse full results for leg
                        full_results = osimresultskey.results["analyses"]["joint_angular_impulse"][foot]
                        
                        # get data
                        implabels = ["net", "windows"]
                        for implabel in implabels:
                            
                            # net impulse
                            if implabel == "net":
                                winlabel = "net"
                                for v, varlabel in enumerate(full_results["net"]["headers"]):
                                    
                                    # skip pelvis and knee beta force
                                    if ("pelvis" in varlabel) or ("knee_angle_beta_force" in varlabel): continue                                    
                                    
                                    # values
                                    net = full_results["net"]["net"][v]
                                    pos = full_results["net"]["pos"][v]
                                    neg = full_results["net"]["neg"][v]

                                    # csv row
                                    csvrow = [subj, group, trial, movement, foot, data_leg_role, implabel, winlabel, varlabel, net, pos, neg]
                                    csvdata.append(csvrow)
                                    
                            elif implabel == "windows":
                                for winlabel in windows:
                                    for v, varlabel in enumerate(full_results["windows"]["headers"]):
                                        
                                        # skip pelvis
                                        if "pelvis" in varlabel: continue
                                        
                                        # values
                                        net = full_results["windows"][winlabel]["net"][v]
                                        pos = full_results["windows"][winlabel]["pos"][v]
                                        neg = full_results["windows"][winlabel]["neg"][v]
    
                                        # csv row
                                        csvrow = [subj, group, trial, movement, foot, data_leg_role, implabel, winlabel, varlabel, net, pos, neg]
                                        csvdata.append(csvrow)                                    
                                    
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                else:
                    print("Dynamic trial: %s" % trial)        
    
    
    # create dataframe
    # appearing in the output
    print("\nCreating dataframe...")
    headers = ["subject", "group", "trial", "movement", "data_leg", "data_leg_role", "period", "window", "variable", "net", "positive", "negative"]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.project.lower() + "_joint_angular_impulse_all_trials.csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    print("\n")                
    
    return None       






'''
-----------------------------------
---- FUNCTIONS: MISCELLANEOUS -----
-----------------------------------
'''


'''
filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
    Filter timeseries data. Raw data can be a list, or an array with rows
    representing time steps and columns as variables. Set cutoff < 0 if
    filtering not required.
''' 
def filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
    
    # return if cutoff < 0 (i.e. filtering not required)
    if cutoff < 0: return data_raw
    
    # filter design
    Wn = sample_rate / 2
    normalised_cutoff = cutoff / Wn
    b, a = signal.butter(butter_order, normalised_cutoff, "lowpass")
    
    # apply filter
    data_filtered = signal.filtfilt(b, a, data_raw, axis = 0)

    return data_filtered


    
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