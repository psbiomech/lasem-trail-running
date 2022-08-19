# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 05:59:34 2022

@author: Prasanna Sritharan
"""



import opensim as osim


modelproc = osim.ModelProcessor("FullBodyModel_SimpleArms_Hamner2010_OpenSim4.3.osim")
modelproc.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
model = modelproc.process()

model.printToXML("FullBodyModel_SimpleArms_Hamner2010_OpenSim4.3_DGF.osim")