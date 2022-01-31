# -*- coding: utf-8 -*-
"""
Run OpenSim pipeline, write input data files

@author: Prasanna Sritharan
"""

import opensim
import pandas as pd
import numpy as np
import pickle as pk
import os
import shutil


'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES




'''
-----------------------------------
---- FUNCTIONS: OPENSIM TOOLS -----
-----------------------------------
'''


'''
opensim_pipeline(meta, user):
    Run OpenSim pipeline: Scale, IK, ID, SO, RRA, CMC.
'''
def opensim_pipeline(meta, user, analyses):
    
    # lower case analyses list
    analyses = [a.casefold() for a in analyses]
    
    # run OpenSim for all valid trials
    for subj in meta:
        for group in meta[subj]["trials"]:
            
            # find the static trial for model scaling
            modelfullpath = ""
            modelfile = ""
            for trial in meta[subj]["trials"][group]:
                if meta[subj]["trials"][group][trial]["usedstatic"]:
                    
                    # load the OsimKey
                    pkpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = meta[subj]["trials"][group][trial]["trial"] + "_osimkey.pkl"
                    with open(os.path.join(pkpath, pkfile),"rb") as fid: 
                        osimkey = pk.load(fid)
                    
                    # run the scale tool if requested
                    if "scale" in analyses:
                        run_opensim_scale(osimkey, user)
                        analyses.remove("scale")
                    
                    # get the full model path
                    modelfile = meta[subj]["trials"][group][trial]["osim"]
                    modelfullpath = os.path.join(pkpath, modelfile)
            
            # find dynamic trials and run requested analyses
            for trial in meta[subj]["trials"][group]:
                if not meta[subj]["trials"][group][trial]["isstatic"]:

                    # load the OsimKey
                    pkpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = meta[subj]["trials"][group][trial]["trial"] + "_osimkey.pkl"
                    with open(os.path.join(pkpath, pkfile),"rb") as fid: 
                        osimkey = pk.load(fid)
                    
                    # copy the model into the trial folder
                    shutil.copy(modelfullpath, pkpath)
                    
                    # run the required analyses
                    for ans in analyses:
                        if not os.path.exists(os.path.join(pkpath, ans)): os.makedirs(os.path.join(pkpath, ans))
                        if ans == "ik":
                            run_opensim_ik(osimkey, user)
                        elif ans == "id":
                            run_opensim_id(osimkey, user)
                        elif ans == "so":
                            run_opensim_so(osimkey, user)
                        elif ans == "rra":
                            #run_opensim_rra(osimkey, user)
                            pass
                        elif ans == "cmc":
                            #run_opensim_cmc(osimkey, user)
                            pass
   
                            
    return None



'''
run_opensim_scale(osimkey, user):    
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run 
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_scale(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    model = osimkey.subject
    trial = osimkey.trial

    print("\nCreating scaled model: %s" % model)
    print("------------------------------------------------")
    
    # create an ScaleTool from a generic setup file
    print("Create new ScaleTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupscale
    tool = opensim.ScaleTool(os.path.join(refsetuppath, refsetupfile))
    tool.setPathToSubject("")
    
    # set subject mass
    tool.setSubjectMass(osimkey.mass)
    
    # set subject name
    tool.setName(osimkey.subject)
    
    
    # ******************************
    # GENERIC MODEL MAKER
    
    print("Initialising GenericModelMaker...")
    
    # set the model file name
    refmodelpath = user.refmodelpath
    refmodelfile = user.refmodelfile
    modelmaker = tool.getGenericModelMaker()
    modelmaker.setModelFileName(os.path.join(refmodelpath, refmodelfile))
    

    # ******************************
    # MODEL SCALER
    
    print("Setting up ModelScaler...")
    
    # set the static TRC file to be used for scaling
    modelscaler = tool.getModelScaler()
    modelscaler.setMarkerFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # set time window
    twindow = opensim.ArrayDouble(0, 2)
    twindow.set(0, 0.50)
    twindow.set(1, 0.55)
    modelscaler.setTimeRange(twindow)
    
    # set output model file name
    modelscaler.setOutputModelFileName(os.path.join(fpath, model + ".osim"))
    
    
    # ******************************
    # MARKER PLACER   

    print("Setting up MarkerPlacer...")

    # set the static TRC file to be used for scaling
    markerplacer = tool.getMarkerPlacer()
    markerplacer.setMarkerFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # set time window
    markerplacer.setTimeRange(twindow)
    
    # set static trial output motion file
    markerplacer.setOutputMotionFileName(os.path.join(fpath, trial + "_static_ik.mot"))
    
    # set output model file
    markerplacer.setOutputModelFileName(os.path.join(fpath, model + ".osim"))    
     
    # set output marker file
    markerplacer.setOutputMarkerFileName(os.path.join(fpath, trial + "_markers.xml"))   


    # ******************************
    # RUN TOOL
    
    print("Running the ScaleTool...")
    
    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + "_Setup_Scale.xml"))

    # run the tool
    try:
        tool.run()
        print("Done.")
    except:
        print("---> ERROR: Scale failed. Skipping Scale for %s." % trial)
    finally:
        print("------------------------------------------------\n")
    
    # ******************************    
    
    return None



'''
run_opensim_ik(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_ik(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    print("Performing IK on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an IK Tool from a generic setup file
    print("Create new IKTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupik
    tool = opensim.InverseKinematicsTool(os.path.join(refsetuppath, refsetupfile))
    
    # load the model
    print("Loading the model: %s..." % modelfile)
    model = opensim.Model(os.path.join(fpath, modelfile))
    model.initSystem()
    
    # set the model in the tool
    tool.setModel(model)

    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][-1])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setStartTime(t0)
    tool.setEndTime(t1)

    # set input TRC file
    print("Setting the marker data (TRC) file...")
    tool.setMarkerDataFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # output MOT file location
    print("Setting output file name...")
    motfilepath = os.path.join(fpath, user.ikcode)
    if not os.path.isdir(motfilepath): os.makedirs(motfilepath)
    tool.setOutputMotionFileName(os.path.join(motfilepath, trial + "_ik.mot"))  


    # ******************************
    # RUN TOOL

    print("Running the IKTool...")

    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + '_Setup_IK.xml'))
        
    # run the tool
    try:
        tool.run()
        print("Done.")
    except:
        print("---> ERROR: IK failed. Skipping IK for %s." % trial)
    finally:
        print("------------------------------------------------\n")
    
    # ******************************
    
    return None



'''
run_opensim_id(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_id(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    print("Performing ID on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an ID Tool from a generic setup file
    print("Create new IDTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupid
    tool = opensim.InverseDynamicsTool(os.path.join(refsetuppath, refsetupfile))

    # load the model
    print("Loading the model: %s..." % modelfile)
    model = opensim.Model(os.path.join(fpath, modelfile))
    model.initSystem()
    
    # set the model in the tool
    tool.setModel(model)   
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][-1])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setStartTime(t0)
    tool.setEndTime(t1)

    # set input directory and coordinates data file
    print("Setting coordinates data file...")
    tool.setInputsDir(fpath)
    tool.setCoordinatesFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(6.0)
    
    # set output directories and generalised forces storage file (note:
    # InverseDynamicsTool XML parser does not seem to like full paths
    # for the OutputGenForceFileName tag, so need to set results dir)
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.idcode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    tool.setOutputGenForceFileName(trial + "_id.sto")
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  

    # set the external loads file name in the inverse dynamics tool
    tool.setExternalLoadsFileName(extloadsfile)
    
    
    # ******************************
    # RUN TOOL 
    
    print("Running the IDTool...")

    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + '_Setup_ID.xml'))
    
    # run the tool
    try:
        tool.run()
        print("Done.")
    except:
        print("---> ERROR: ID failed. Skipping ID for %s." % trial)
    finally:
        print("------------------------------------------------\n")

    # ******************************
        
    return None    
    


'''
run_opensim_so(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_so(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    print("Performing SO on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an generic AnalyzeTool from the template setup file (set the flag
    # aLoadModelAndInput = false as we are only configuring the setup file at 
    # this stage)
    print("Create new AnalyzeTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupso
    tool = opensim.AnalyzeTool(os.path.join(refsetuppath, refsetupfile), False)
    tool.setName(trial)
    
    # set the model in the tool
    print("Loading the model: %s..." % modelfile)
    tool.setModelFilename(os.path.join(fpath, modelfile))   
 
    # append reserve actuators (need to create a new ArrayStr with one element, 
    # the file name, and then pass the ArrayStr to setForceFiles)
    print("Appending reserve actuators...")
    tool.setReplaceForceSet(False)
    fsvec = opensim.ArrayStr()
    refrefreserveactuators = user.refreserveactuators
    fsvec.append(os.path.join(refsetuppath, refrefreserveactuators))
    tool.setForceSetFiles(fsvec)
        
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][-1])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)

    # set coordinates data file
    print("Setting coordinates data file...")
    tool.setCoordinatesFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(6.0)
    
    # set output directories and generalised forces storage file (note:
    # InverseDynamicsTool XML parser does not seem to like full paths
    # for the OutputGenForceFileName tag, so need to set results dir)
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.socode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  

    # set the external loads file name in the inverse dynamics tool
    tool.setExternalLoadsFileName(extloadsfile)
    
    # create an SO analysis
    print("Append new SO Analysis to the AnalysisSet...")
    so = opensim.StaticOptimization()
    so.setName("so")
    so.setStartTime(t0)
    so.setEndTime(t1)
    
    # add the StaticOptimization analysis to the tool AnalysisSet
    analyses = tool.getAnalysisSet()
    analyses.insert(0, so)
    
    
    # ******************************
    # RUN TOOL 
    
    print("Running the AnalysisTool (SO)...")

    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + '_Setup_SO.xml')
    tool.printToXML(customsetupfile)
    
    # run the tool (need to load the setup again into a new AnalyzeTool)
    try:
        tool2 = opensim.AnalyzeTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: SO failed. Skipping SO for %s." % trial)
    finally:
        print("------------------------------------------------\n")

    # ******************************
        
    return None    
    

        

'''
-----------------------------------
----- FUNCTIONS: OPENSIM DATA -----
-----------------------------------
'''


'''
write_ground_forces_mot_file(osimkey):
    Write .mot file for ground forces (external loads).
'''
def write_ground_forces_mot_file(osimkey):

    # output dataframe info
    ns = len(osimkey.forces["time"])
    nc = 19    
    t0 = osimkey.forces["time"][0]
    tf = osimkey.forces["time"][-1]

    # write headers
    fname = osimkey.trial + "_grf.mot"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("%s\n" % fname)
        f.write("nRows=%d\n" % ns)
        f.write("nColumns=%s\n" % nc)
        f.write("\n")
        f.write("name %s\n" % fname)
        f.write("datacolumns %d\n" % nc)
        f.write("datarows %d\n" % ns)
        f.write("range %f %f\n" % (t0, tf))
        f.write("endheader\n")

    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.forces["time"]
    datamat[:,1:4] = osimkey.forces["data"]["right"]["F"]
    datamat[:,4:7] = osimkey.forces["data"]["right"]["cop"]
    datamat[:,7:10] = osimkey.forces["data"]["left"]["F"]
    datamat[:,10:13] = osimkey.forces["data"]["left"]["cop"]
    datamat[:,13:16] = osimkey.forces["data"]["right"]["T"]
    datamat[:,16:19] = osimkey.forces["data"]["left"]["T"]    
        
    # convert to dataframe
    headers = ["time", "grf_right_vx", "grf_right_vy", "grf_right_vz", "grf_right_px", "grf_right_py", "grf_right_pz", "grf_left_vx", "grf_left_vy", "grf_left_vz", "grf_left_px", "grf_left_py", "grf_left_pz", "grf_right_tx", "grf_right_ty", "grf_right_tz", "grf_left_tx", "grf_left_ty", "grf_left_tz"]
    data = pd.DataFrame(datamat, columns=headers)
        
    # write table
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=True, index=False)
    
    return data



'''
write_marker_trajctory_trc_file(osimkey):
    Write .trc file for marker trajectories.
'''
def write_marker_trajctory_trc_file(osimkey):
    
    # output dataframe info
    ns = len(osimkey.markers["time"])
    nm = len(osimkey.markers) - 5
    nc = 2 + (nm * 3)
    rate = osimkey.markers["rate"]
    units = osimkey.markers["units"]

    # remove non-marker dict keys
    markernames0 = list(osimkey.markers.keys())
    markernames0.remove("rate")
    markernames0.remove("units")
    markernames0.remove("offset")
    markernames0.remove("frames")
    markernames0.remove("time")

    # build marker headers
    markernames = ""
    markernames = markernames.join(["Frame#\t","Time\t"] + list(map(lambda x: x + "\t\t\t", markernames0)))
    dirnums = list(range(1,len(markernames0) + 1))
    dirnames = ""
    dirnames = dirnames.join(["\t"] + list(map(lambda n: "\tX%d\tY%d\tZ%d" % (n, n, n), dirnums)))

    # write headers
    fname = osimkey.trial + "_markers.trc"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("PathFileType\t4\t(X/Y/Z)\t%s\n" % fname)
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        f.write("%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (rate, rate, ns, nm, "mm", rate, 1, ns))
        f.write("%s\n" % markernames)
        f.write("%s\n" % dirnames)
        
    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.markers["frames"]
    datamat[:,1] = osimkey.markers["time"]
    n = 2
    for mkr in markernames0:
        mkrdata = osimkey.markers[mkr]
        if units.casefold() == "m": mkrdata = mkrdata * 1000
        datamat[:,n:n+3] = mkrdata
        n = n + 3
    
    # convert to dataframe
    data = pd.DataFrame(datamat)
    data[0] = data[0].astype(int)
    
    # write table, no headers
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=False, index=False)
    
    return data