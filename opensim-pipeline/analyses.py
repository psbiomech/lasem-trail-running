# -*- coding: utf-8 -*-
"""
OpenSim post-hoc analyses

@author: Prasanna Sritharan
"""



import os
import numpy as np
import pickle as pk
import pandas as pd
from scipy import interpolate
from scipy import integrate



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
analyses_batch_process(meta, user):
    Batch process post-hoc analyses.
'''
def analyses_batch_process(meta, user, analyses):
    
    # extract OpenSim data
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

                    # Load OsimResultsKey
                    pklpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        osimresultskey = pk.load(fid0)
                            
                    # Output analysis dict
                    analyses = {}
                    analyses["subj"] = subj
                    analyses["group"] = group
                    analyses["trial"] = trial
                    analyses["outpath"] = pklpath
                                  
                    print("Dynamic trial: %s" % trial)
                    
                    
                    
                    # ******************************
                    # ANALYSES
    
                    for ans in analyses:
                        
                        # Joint angular impulse
                        if ans == "jai":
                            print("---> Joint angular impulse")
                            analyses["joint_angular_impulse"] = calculate_joint_angular_impulse(osimresultskey, user)     
                            
                        # Joint angular power
                        elif ans == "jap":
                            print("---> Joint angular power")
                            analyses["joint_angular_power"] = calculate_joint_angular_power(osimresultskey, user)
                            
                        # Joint angular work
                        elif ans == "jaw":
                            print("---> Joint angular work")
                            analyses["joint_angular_work"] = calculate_joint_angular_work(osimresultskey, user)                            
                            
        
    
                    #                
                    # ******************************
     
                    # pickle analysis dict
                    with open(os.path.join(pklpath, trial + "_analyses_results.pkl"), "wb") as fid1:
                        pk.dump(osimresultskey, fid1)   

                except:
                    print("*** FAILED ***")
                    failedfiles.append(trial)                    
                        
                          
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
calculate_joint_angular_power(osimresultskey, user):
    Calculate the joint angular power for all available joints using joint angles
    from inverse kinematics and joint moments from inverse dynamics.
'''
def calculate_joint_angular_power(osimresultskey, user):
    
    # Get event times
    events = osimresultskey.events["time"]
    
    # Get ID data
    Tdata = osimresultskey.results["split"][user.idcode].copy()
    
    # Get IK data
    Qdata = osimresultskey.results["split"][user.ikcode].copy()
    
    # Calculate joint angular power on each leg
    jap = {}
    for leg in user.leg:
        
        # Leg joint moments (Nm), remove pelvis forces
        timeT = Tdata[leg]["data"][:, 0]
        Torig = Tdata[leg]["data"][:, 1:]
        headersT = Tdata[leg]["headers"][1:]
        delidxT = headersT.index("pelvis_tx_force")
        del headersT[delidxT:delidxT+3]
        Torig = np.delete(Torig, range(delidxT, delidxT + 3), axis = 1)       
               
        # Leg joint angles (rad), remove pelvis translations
        timeQ = Qdata[leg]["data"][:, 0]
        Q = np.deg2rad(Qdata[leg]["data"][:, 1:])
        headersQ = Qdata[leg]["headers"][1:]
        delidxQ = headersQ.index("pelvis_tx")
        del headersQ[delidxQ:delidxQ+3]
        Q = np.delete(Q, range(delidxQ, delidxQ + 3), axis = 1)
        
        # Calculate joint velocities (rad/s)
        Qdot = np.gradient(Q, timeQ, axis = 0)
        
        # Adjust moments using interpolation so that time vector matches
        # the velocity data. OpenSim can cause small shifts in the time
        # vectors that need to be corrected.
        Tinterpfun = interpolate.interp1d(timeT, Torig, kind = "cubic", axis = 0, fill_value = "extrapolate")
        T = Tinterpfun(timeQ)
        
        # Calculate joint angular power for each coordinate
        P = np.zeros(np.shape(Q))
        for idxQ, qcoord in enumerate(headersQ):
            
            # Find moment column index of current coordinate as moments and
            # velocities may be out of order
            idxT = [1 if qcoord in head else 0 for i, head in enumerate(headersT)].index(1)
            
            # Get the data for the current coordinate
            qdot = Qdot[:, idxQ]
            t = T[:, idxT]

            # Calculate joint angular power: P = Tw (Nm/s)
            P[:, idxQ] = np.multiply(t, qdot)
            
        # Add time column
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
    
    # Calculate the joint angular work
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



    
