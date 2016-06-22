
# coding: utf-8

# In[6]:

import argparse
import os
import numpy as np
import glob

def valid_file(parser,arg):
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file does not exist")
    else:
        return arg

parser = argparse.ArgumentParser(description="Input an mseID and the freeview command with relevant information will output")
parser.add_argument("-m","--mse",
                    dest="mse",
                    type = str,
                    default=[],
                    metavar="mse####",
                    action="append",
                    help="import one or more mseIDs using the format mse####")
args = parser.parse_args()

if args.mse == []:
    raise parser.error("No mseIDs inputted")
else:
    subjects = args.mse

for num in range(len(subjects)):
    for les_file in glob.glob("/data/henry7/PBR/subjects/%s/lesions_manual/*/alignment_lesions.nii.gz" % subjects[num]):
        les_file
    for seg_file in glob.glob("/data/henry7/PBR/subjects/%s/masks/*/segmentation.nii.gz" % subjects[num]):
        seg_file
    for gm_file in glob.glob("/data/henry6/PBR/surfaces/*%s*/mri/ribbon.mgz" % subjects[num]):
        gm_file
    print; print "Freeview command for %s" % subjects[num]
    print; print "freeview",gm_file,seg_file,les_file; print


# In[ ]:



