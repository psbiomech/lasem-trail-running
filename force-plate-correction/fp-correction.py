# -*- coding: utf-8 -*-
"""
Calculate and apply offset to no-zeroed FP data

@author: Prasanna Sritharan
"""

import os
import glob
import datetime
import csv
import numpy as np
import pyc3dserver as c3d



print("\n\n\n")
print("----------------------------------------")
print("TRAIL: FORCE PLATE ZERO CORRECTION")
print("----------------------------------------")

# start time stamp
ts0 = datetime.datetime.now();
print("Start: %s" % ts0)

print("----------------------------------------")
print("\n")



# Load C3DServer
itf = c3d.c3dserver()


print("\n\nCalculating offsets using all dynamic C3D files...")

# Get all trials for subject and estimate offset for each channel
filepath = r"C:\Users\Owner\Documents\data\TRAIL Test Data\fp-test"
c3dlist = glob.glob(os.path.join(filepath, "*.c3d"))
offset = {}
channel_prefixes = ["Force", "Moment"]
for filename in c3dlist:

    # Skip static trials
    if "static" in filename.casefold(): continue    

    # Get the raw FP analog streams
    c3d.open_c3d(itf, os.path.join(filename))
    fp = c3d.get_dict_analogs(itf, excl_forces=False)
    
    # Find average force plate offset from all dynamic trials
    sample_window = 100
    for label in fp["LABELS"]:
        if any([prefix in label for prefix in channel_prefixes]):
            data0 = c3d.get_analog_data_unscaled(itf, label)
            if label in offset:
                offset[label] = np.append(offset[label], np.mean(data0[-sample_window]))
            else:
                offset[label] = np.array(np.mean(data0[-sample_window]))
       
    # Close the C3D file
    c3d.close_c3d(itf)
      
        
      
print("Calculating mean offsets for each analog stream...")  

# Calculate mean offsets for that subject
mean_offset = {}
for label in offset:
    mean_offset[label] = np.mean(offset[label])
    


print("Applying mean offsets to all trials, including static trials...")  

# Apply offsets to all trials for that subject, including static trials, and
# save the corrected file
outdir = os.path.join(filepath, "offset")
if not(os.path.exists(outdir)): os.makedirs(outdir)
for filename in c3dlist:
    c3dfile = c3d.open_c3d(itf, os.path.join(filename))
    for label in mean_offset:
        data0 = c3d.get_analog_data_unscaled(itf, label)
        data1 = data0 - mean_offset[label]
        c3d.set_analog_data(itf, label, data1)
    c3d.save_c3d(itf, os.path.join(filepath, outdir, os.path.splitext(os.path.basename(filename))[0] + "_offset.c3d"))
    c3d.close_c3d(itf)



print("Writing offsets to CSV file...")  

# Write offsets to file
with open(os.path.join(filepath, outdir, "TRAIL_324_offsets.csv"), "w", newline='') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in mean_offset.items():
       writer.writerow([key, value])



    
print("\n")
print("----------------------------------------")

# end time stamp
ts1 = datetime.datetime.now();
print("End: %s" % ts1)
datetime.datetime.now()

print("----------------------------------------")
print("\n")
