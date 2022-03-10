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

    # clear OpenSim log file
    #logfile0 = os.path.join(user.logfilepath, user.logfile)    
    #if os.path.isfile(logfile0): os.remove(logfile0)
    
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
                    
                    # copy the OpenSim log file into the local folder
                    #logfile1 = os.path.join(pkpath, meta[subj]["trials"][group][trial]["trial"] + ".log")
                    #shutil.copyfile(logfile0, logfile1)                    
            
            # find dynamic trials and run requested analyses
            for trial in meta[subj]["trials"][group]:
                
                
                # ****** FOR TESTING ONLY ******
                # import re
                # trialre = re.compile("TRAIL_071_EP_08")
                # if not trialre.match(trial):
                #     print("%s ---> SKIP" % trial)
                #     continue
                # ******************************
                
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
                            run_opensim_rra(osimkey, user)
                        elif ans == "cmc":
                            run_opensim_cmc(osimkey, user)
                            
                    # copy the OpenSim log file into the local folder
                    #logfile1 = os.path.join(pkpath, meta[subj]["trials"][group][trial]["trial"] + ".log")
                    #shutil.copyfile(logfile0, logfile1)  
   
                            
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
    subject = osimkey.subject
    model = osimkey.model
    trial = osimkey.trial

    print("\nCreating scaled model: %s" % subject)
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
    modelscaler.setOutputModelFileName(os.path.join(fpath, model))
    
    
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
    markerplacer.setOutputModelFileName(os.path.join(fpath, model))    
     
    # set output marker file
    markerplacer.setOutputMarkerFileName(os.path.join(fpath, trial + "_markers.xml"))   


    # ******************************
    # RUN TOOL
    
    print("Running the ScaleTool...")
    
    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + "_Setup_Scale.xml"))

    # run the tool, and scale the model if required
    try:
        
        # run the tool
        tool.run()
        
        # scale the model FoM if required
        sf_fom = user.fom_scalefactor
        if (type(sf_fom) is dict) or (sf_fom >= 0):
            print("---> Scaling muscle FoM in model...")
            shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_FoM.osim"))
            model0 = opensim.Model(os.path.join(fpath, model))
            refmodel = opensim.Model(os.path.join(refmodelpath, refmodelfile))
            model1 = update_osim_fom(model0, sf_fom, refmodel)
            model1.printToXML(os.path.join(fpath, model))
            
        # scale the model FoM if required
        sf_lst = user.lst_scalefactor
        if (type(sf_lst) is dict) or (sf_lst > 0):
            print("---> Scaling muscle LsT in model...")
            shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_LsT.osim"))
            model0 = opensim.Model(os.path.join(fpath, model))
            model1 = update_osim_lst(model0, sf_lst)
            model1.printToXML(os.path.join(fpath, model))            
                    
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

    # set the initial and final times (limit to between first and last event,
    # but allow extra 0.05 sec at start for CMC)
    t0 = float(osimkey.events["time"][0]) + user.cmc_start_time_offset
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
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
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

    # set the external loads file name
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
  
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)

    # set coordinates data file
    print("Setting coordinates data file...")
    tool.setCoordinatesFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(6.0)
    
    # set output directory
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

    # set the external loads file name
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
    # PREPARE RESERVE ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
        
    # get the pelvis COM location from the model
    model = opensim.Model(os.path.join(fpath, modelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter() 
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refrefreserveactuators = user.refreserveactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, refrefreserveactuators))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write updated actuator set to file    
    forcesetfile = trial + "_Reserve_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
      
    # append reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    print("Appending reserve actuators...")
    fsvec = opensim.ArrayStr()    
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(False)   

    # # append default reference actuators, with pelvis residuals acting at 
    # # ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refreserveactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_Reserve_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(False)     

    
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
run_opensim_rra(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_rra(osimkey, user):

    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial
    
    print("Performing RRA on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an RRA Tool from a generic setup file
    # (note: set loadModel=false for manual load later)
    print("Create new RRATool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetuprra
    tool = opensim.RRATool(os.path.join(refsetuppath, refsetupfile), False)      
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0]) + user.cmc_start_time_offset
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]]) + user.rra_end_time_offset
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)

    # set flag to adjust model COM and set body to torso
    tool.setAdjustCOMToReduceResiduals(True)
    tool.setAdjustedCOMBody("torso")

    # set desired kinematics file name (original IK results)
    tool.setDesiredKinematicsFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    
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
    
    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    # set output directory
    print("Setting output folder name...")
    stofilepath = os.path.join(fpath, user.rracode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)


    # ******************************
    # PREPARE RRA ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
    
    # get the pelvis COM location from the model
    model = opensim.Model(os.path.join(fpath, modelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter()
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refforcesetfile = user.refrraactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write actuator set to file
    forcesetfile = trial + "_RRA_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
    
    # replace reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    fsvec = opensim.ArrayStr();
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(True)
    
    # # replace actuators with default reference actuators, with pelvis residuals
    # # acting at ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refrraactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_RRA_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))    
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(True)    
    
    
    # ******************************
    # PREPARE RRA TASKS
 
    # load RRA tasks set
    rratasksfile = user.refrratasks
    rrataskset = opensim.CMC_TaskSet(os.path.join(refsetuppath, rratasksfile))
    
    # updates here...
    # e.g. controller weights and gains
    
    # print tasks to trial folder
    rrataskset.printToXML(os.path.join(fpath, trial + "_RRA_Tasks.xml"))
    
    # set RRA tasks file in tool
    tool.setTaskSetFileName(os.path.join(fpath, trial + "_RRA_Tasks.xml"))
       
    
    # ******************************
    # RUN TOOL
    
    print("Running the RRATool, %d iterations...\n---> (this may take a while)" %  user.rraiter)
    
    try:

        # run the tool require number of iterations
        rraadjustedmodel = []
        for i in range(user.rraiter):
            
            print("---> Iteration: %d" % (i + 1))
            
            # set tool name based on current iteration
            rraiter = "RRA_" + str(i + 1)
            rraname = trial + "_" + rraiter
            tool.setName(rraname)
            
            # get the model, ensure model MTP joints are locked (only needs to
            # be done on the first iteration
            if i==0:
                rramodelfile = os.path.join(fpath, modelfile)
                #lockMTPJoints(rramodelfile)
            else:
                rramodelfile = rraadjustedmodel 
                
            # set model file name
            tool.setModelFilename(rramodelfile)
            
            # set new adjusted model file name
            rraadjustedmodel = os.path.join(fpath, rraname + "_AdjustedModel.osim")
            tool.setOutputModelFileName(rraadjustedmodel)
            
            # save the current settings in a setup file
            rrasetupfile = os.path.join(fpath, trial + "_Setup_" + rraiter + ".xml")
            tool.printToXML(rrasetupfile)
            
            # load the current setup file, run the current RRA tool
            rratool2 = opensim.RRATool(rrasetupfile);
            rratool2.run();            
        
        print("Done.")
        
    except:
        print("---> ERROR: RRA failed. Skipping RRA for %s." % trial)
    finally:
        print("------------------------------------------------\n")
    
    # ******************************
        
    return None    
        


'''
run_opensim_cmc(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_cmc(osimkey, user):

    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial
    
    print("Performing CMC on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create a CMC Tool from a generic setup file
    # (note: set loadModel=false for manual load later)
    print("Create new CMCTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupcmc
    tool = opensim.CMCTool(os.path.join(refsetuppath, refsetupfile), False)      
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0]) + user.cmc_start_time_offset
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    
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
    
    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    # set output directory
    print("Setting output folder name...")
    stofilepath = os.path.join(fpath, user.cmccode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)

    # use fast target (should normally be set to true, but set to false if 
    # having difficulty solving, e.g. pathological movement patterns or dynamic
    # movement patterns such as fast running)
    tool.setUseFastTarget(user.use_fast_target)


    # ******************************
    # PREPARE CMC ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
    
    # get the pelvis COM location from the desired model
    # (the pelvis COM should be the same for both models anyway)
    if user.use_rra_model:
        rramodelfile = trial + "_RRA_" + str(user.rraiter) + "_AdjustedModel.osim"
    else:
        rramodelfile = modelfile   
    model = opensim.Model(os.path.join(fpath, rramodelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter()
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refforcesetfile = user.refcmcactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write actuator set to file
    forcesetfile = trial + "_CMC_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
    
    # replace reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    fsvec = opensim.ArrayStr();
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(False)
    
    # # append default reference actuators, with pelvis residuals acting at 
    # # ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refcmcactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_CMC_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))      
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(False)    
    
    
    # ******************************
    # PREPARE CMC TASKS
    # (allows for trial-specific controller weights and gains to be applied)
 
    # load CMC tasks set
    cmctasksfile = user.refcmctasks
    cmctaskset = opensim.CMC_TaskSet(os.path.join(refsetuppath, cmctasksfile))
    
    # updates here... (TBD)
    # e.g. controller weights and gains
    
    # print tasks to trial folder
    cmctaskset.printToXML(os.path.join(fpath, trial + "_CMC_Tasks.xml"))
    
    # set CMC tasks file in tool
    tool.setTaskSetFileName(os.path.join(fpath, trial + "_CMC_Tasks.xml"))
       

    # ******************************
    # PREPARE CMC CONTROL CONSTRAINTS
    # (allows for trial-specific controller constraints to be applied)

    # load reference CMC control constraints
    cmccontrolsfile = user.refcmccontrolconstraints
    cmccontrolset = opensim.ControlSet(os.path.join(refsetuppath, cmccontrolsfile))
    
    # updates here... (TBD)
    # e.g. actuator gains
    
    # print control constraints to trial folder
    cmccontrolset.printToXML(os.path.join(fpath, trial + '_CMC_ControlConstraints.xml'))
    
    # set CMC control constraints file in tool
    tool.setConstraintsFileName(os.path.join(fpath, trial + '_CMC_ControlConstraints.xml'))
    

    # ******************************
    # SET THE MODEL AND KINEMATICS
    # (use baseline or RRA adjusted model and kinematics)

    # set the desired model
    if user.use_rra_model:
        actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_AdjustedModel.osim"
    else:
        actualmodelfile = modelfile
    tool.setModelFilename(os.path.join(fpath, actualmodelfile))
        
    # set desired kinematics file and filter frequency
    if user.use_rra_model:
        kinfile = os.path.join(fpath, user.rracode, trial + "_RRA_" + str(user.rraiter) + "_Kinematics_q.sto")
        filtfreq = -1
    else:
        kinfile = os.path.join(fpath, user.ikcode, trial + "_ik.mot")    
        filtfreq = 6.0
    tool.setDesiredKinematicsFileName(kinfile)
    tool.setLowpassCutoffFrequency(filtfreq)
   
    
    # ******************************
    # RUN TOOL
    
    print("Running the CMCTool...\n---> (this may take a while)")
 
    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + '_Setup_CMC.xml')
    tool.printToXML(customsetupfile)   
 
    # run the tool (need to load the setup again into a new CMCTool and set the
    # use_fast_target flag)
    try:
        tool2 = opensim.CMCTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: CMC failed. Skipping CMC for %s." % trial)
    finally:
        print("------------------------------------------------\n")
    
    # ******************************
        
    return None    
    
    


'''
-----------------------------------
---- FUNCTIONS: MISCELLANEOUS -----
-----------------------------------
'''



'''
update_osim_fom(modelfullpath, scalefactor, refmodelpath):
    Update all OpenSim model FoM by a fixed scale factor, or using the scaling
    law described by Handsfield et al. 2013 (scalefactor["all"] = -1), and/or 
    apply a custom scale factor to selected muscles. To ignore scale factors,
    set scalefactor["all"] = 0.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", -1): apply Handsfield scaling law
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0        
'''
def update_osim_fom(model, scalefactor, refmodel):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale by a fixed scale factor, or use Handsfield scaling law
    sf = 1.0
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                currmuscle = allmuscles.get(m)
                currmuscle.setMaxIsometricForce(sf * currmuscle.getMaxIsometricForce())
        elif scalefactor["all"] == -1:
            # Handsfield scaling law: TBD
            return model
        elif scalefactor["all"] == 0:
            return model
                        
    # custom scale selected variables, overwrites muscles scaled by "all" key
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                sfm = scalefactor[sfname] / sf
                currmuscle.setMaxIsometricForce(sfm * currmuscle.getMaxIsometricForce())

        
    return model
            


'''
update_osim_lst(modelfullpath, scalefactor):
    Update all OpenSim model LsT by a fixed scale factor, and/or apply a
    custom scale factor to selected muscles. To ignore scale factors, set
    scalefactor["all"] = 0.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0
'''
def update_osim_lst(model, scalefactor):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale all by a fixed scale factor
    sf = 1.0
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                currmuscle = allmuscles.get(m)
                currmuscle.set_tendon_slack_length(sf * currmuscle.get_tendon_slack_length())
        elif scalefactor["all"] == 0:
            return model

    # custom scale selected variables
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                sfm = scalefactor[sfname] / sf
                currmuscle.set_tendon_slack_length(sfm * currmuscle.get_tendon_slack_length())
                    
    return model
            


'''
write_ground_forces_mot_file(osimkey):
    Write .mot file for ground forces (external loads).
'''
def write_ground_forces_mot_file(osimkey):

    # output dataframe info
    ns = len(osimkey.forces["time"])
    nc = 19    

    # write headers
    fname = osimkey.trial + "_grf.mot"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("%s\n" % fname)
        f.write("version=1\n")
        f.write("nRows=%d\n" % ns)
        f.write("nColumns=%s\n" % nc)
        f.write("inDegrees=yes\n")
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
    headers = ["time", "ground_force_vx", "ground_force_vy", "ground_force_vz", "ground_force_px", "ground_force_py", "ground_force_pz", "1_ground_force_vx", "1_ground_force_vy", "1_ground_force_vz", "1_ground_force_px", "1_ground_force_py", "1_ground_force_pz", "ground_torque_x", "ground_torque_y", "ground_torque_z", "1_ground_torque_x", "1_ground_torque_y", "1_ground_torque_z"]
    data = pd.DataFrame(datamat, columns=headers)
        
    # write table
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=True, index=False, float_format="%20.10f")
    
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
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=False, index=False, float_format="%20.10f")
    
    return data