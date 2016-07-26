
# coding: utf-8

# # The code

# In[8]:

## Command imports
import argparse
import os
import numpy as np

def valid_file(parser,arg):
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file does not exist")
    else:
        return arg

parser = argparse.ArgumentParser(description='Input a mseID or a list (using .txt files) of mseIDs, and if the exam is segmented, a list of lesions and some information for every scan will output in a csv file')
group = parser.add_mutually_exclusive_group()
group.add_argument("-m","--mse", 
                    dest="mseID",
                    type = str,
                    default=[],
                    metavar="mse####",
                    action="append",
                    help="import one mseID using the format 'mse####' (do not add 0 if mse has only 3 digits). repeat if multiple mseIDs are necessary")
group.add_argument("-f", "--textfile", 
                    dest="filename", 
                    default ="no_file_chosen",
                    type=str,
                    help="import mseIDs line-by-line from a txt file")
parser.add_argument("-d", "--decimals",
                    dest="dec",
                    type=int,
                    default = 3,
                    help="indicate how many digits after the decimal point to show (default is 3)")

args = parser.parse_args()

decval = args.dec #number of decimals to round to

if args.mseID == []:
    if args.filename == "no_file_chosen":
        raise parser.error("No mseIDs inputted")
    else:
        msetxt = valid_file(parser,args.filename)
        subjects_txt = np.genfromtxt(msetxt, dtype=str)
        subjects = subjects_txt
elif args.mseID != []:
    subjects_man = args.mseID
    subjects = subjects_man


import nibabel as nib
import glob
import math
import pandas as pd
from scipy.ndimage import label
from numpy.linalg import inv

## Function that gets coordinates from labels

def get_label_coord(labels,num_labels):
    all_label_coords = []
    if num_labels >= 1:
        for count in range(1, num_labels+1):
            cur_label_coords = []
            x,y,z = np.nonzero(labels==count)
            for count2 in range(len(x)):
                cur_coord = [x[count2],y[count2],z[count2]]
                cur_label_coords.append(cur_coord)
            all_label_coords.append(cur_label_coords)
    else:
        x,y,z = np.nonzero(labels==1)
        for count in range(len(x)):
            all_label_coords.append([x[count],y[count],z[count]])
    return all_label_coords

## Function that converts coordinates to real space

def get_rs_coord(coordinates,affine):
    all_rs_coords = []
    for count in range(len(coordinates)):
        cur_rs_coords = []
        for count2 in range(len(coordinates[count])):
            rs_coord = np.dot(affine, [coordinates[count][count2][0],
                                           coordinates[count][count2][1], 
                                           coordinates[count][count2][2],
                                           1])
            rs_coord_noone = [rs_coord[0],rs_coord[1],rs_coord[2]]
            cur_rs_coords.append(rs_coord_noone)
        all_rs_coords.append(cur_rs_coords)
    return all_rs_coords

## Function that converts coordinates to segmentation file coordinates

def get_seg_coord(coordinates):
    all_seg_coords = []
    for count in range(len(coordinates)):
        cur_seg_coords = []
        for count2 in range(len(coordinates[count])):
            seg_coord = np.dot(inv_seg_affine, [coordinates[count][count2][0],
                                                coordinates[count][count2][1], 
                                                coordinates[count][count2][2],
                                                1])
            seg_coord_noone = [int(seg_coord[0]),int(seg_coord[1]),int(seg_coord[2])]
            cur_seg_coords.append(seg_coord_noone)
        all_seg_coords.append(cur_seg_coords)
    return all_seg_coords
    
## Function that takes average of lesion coordinates

def average_func(coordinates):
    sumx=0;sumy=0;sumz=0
    for count in range(len(coordinates)):
        sumx += coordinates[count][0]
        sumy += coordinates[count][1]
        sumz += coordinates[count][2]
    average_x = sumx / len(coordinates)
    average_y = sumy / len(coordinates)
    average_z = sumz / len(coordinates)
    #print count, averages
    return [average_x, average_y, average_z]

## Euclidean distance determination function 

def dist_det(struc_coord, lesion_coord):
    dist_list = []
    for x in range(len(struc_coord)):
        dist_list.append(math.sqrt((struc_coord[x][0] - lesion_coord[0]) ** 2 +
                                   (struc_coord[x][1] - lesion_coord[1]) ** 2 +
                                   (struc_coord[x][2] - lesion_coord[2]) ** 2))
    return dist_list

## Import mse IDs 

all_sub_results = {"mseID":[],
                   "total number of lesions":[],      
                   "subcortical lesions":[],
                   "juxtacortical lesions":[],
                   "periventricular lesions":[],
                   "infratentorial lesions":[],
                   "lesion":[],
                   "type":[],
                   "center coordinates":[],
                   "volume":[],
                   "distance from midbrain":[],
                   "distance from ventricles":[],
                   "distance from gray matter":[]}
problems = []

## Begin Program

for numsubjects in range(len(subjects)):

    ## Set variables to lesion/segmentation files, freeview cmdline copy and paste, makes sure segmentations exist before proceeding through rest of program
    print "Preparing", subjects[numsubjects];print
    
    les_src = glob.glob("/data/henry6/PBR/surfaces/*%s*/recon_edits/aseg_infra.nii.gz" % subjects[numsubjects])
    if len(les_src) != 0:
        les_file = les_src[0]
        les_label_num = 77
    else:
        les_src = glob.glob("/data/henry7/PBR/subjects/%s/lesions_manual/*/alignment_lesions.nii.gz" % subjects[numsubjects])
        if len(les_src) != 0:
            les_file = les_src[0]
            les_label_num = 1
        else:
            les_src = glob.glob("/data/henry7/PBR/subjects/%s/lesions_manual/*/*lesions.nii.gz" % subjects[numsubjects])
            if len(les_src) != 0:
                les_file = les_src[0]
                les_label_num = 1
            else:
                les_src = glob.glob("/data/henry6/PBR/surfaces/*%s*/mri/aseg.mgz" % subjects[numsubjects])
                if len(les_src) != 0:
                    les_file = les_src[0]
                    les_label_num = 77
                else:
                    print "Could not find lesion segmentation file, skipping %s" % subjects[numsubjects]; print
                    problems.append(subjects[numsubjects])
                    continue
    seg_src = glob.glob("/data/henry7/PBR/subjects/%s/masks/*/segmentation.nii.gz" % subjects[numsubjects])
    if len(seg_src) != 0:
        seg_file = seg_src[0]
    else:
        print "Could not find brain segmentation file, skipping %s" % subjects[numsubjects]; print
        problems.append(subjects[numsubjects])
        continue
    gm_src = glob.glob("/data/henry6/PBR/surfaces/*%s*/mri/ribbon.mgz" % subjects[numsubjects])
    if len(gm_src) != 0:
        gm_file = gm_src[0]
    else:
        print "Could not find gray matter segmentation file, skipping %s" % subjects[numsubjects]; print
        problems.append(subjects[numsubjects])
        continue
        
    print "Freeview link:"; print 
    print "freeview", gm_file, seg_file, les_file; print

    ## Obtain affines 

    # les corresponds to lesion file
    # 
    # seg corresponds to segmentation file (important for midbrain and ventricular structures, plus real space conversion)
    #
    # gm corresponds to ribbon file that gives gray matter coordinates

    les_img = nib.load(les_file)
    les_img.dataobj
    seg_img = nib.load(seg_file)
    seg_img.dataobj
    gm_img = nib.load(gm_file)
    gm_img.dataobj

    les_data = les_img.get_data()
    seg_data = seg_img.get_data()
    gm_data = gm_img.get_data()

    les_affine = les_img.get_affine()
    gm_affine = gm_img.get_affine()
    seg_affine = seg_img.get_affine()
    inv_seg_affine = np.linalg.inv(seg_affine)

    ## Set labels for structure's coordinates

    #lesion labels
    les_labels, n_les_labels = label(les_data==les_label_num)

    #segmentation - brainstem labels
    seg_brainstem_labels, n_seg_brainstem_labels = label(seg_data==[16])
    seg_lcerebellumcortex_labels, n_seg_lcerebellumcortex_labels = label(seg_data==[8])
    seg_rcerebellumcortex_labels, n_seg_rcerebellumcortex_labels = label(seg_data==[47])
    seg_lcerebellumwm_labels, n_seg_lcerebellumwm_labels = label(seg_data==[7])
    seg_rcerebellumwm_labels, n_seg_rcerebellumwm_labels = label(seg_data==[46])

    #segmentation - lateral ventricle labels
    seg_llv_labels, n_seg_llv_labels = label(seg_data==[4])
    seg_rlv_labels, n_seg_rlv_labels = label(seg_data==[43])

    #gm labels
    gm_lh_labels, n_gm_lh_labels = label(gm_data==3)
    gm_rh_labels, n_gm_rh_labels = label(gm_data==42)

    ## Generate midbrain coordinates into a variable

    bs = get_label_coord(seg_brainstem_labels,n_seg_brainstem_labels)
    lcc = get_label_coord(seg_lcerebellumcortex_labels,n_seg_lcerebellumcortex_labels)
    lcw = get_label_coord(seg_lcerebellumwm_labels,n_seg_lcerebellumwm_labels)
    rcc = get_label_coord(seg_rcerebellumcortex_labels,n_seg_rcerebellumcortex_labels)
    rcw = get_label_coord(seg_rcerebellumwm_labels,n_seg_rcerebellumwm_labels)

    midbrain_tmp = []
    midbrain = []

    if len(bs) != 1:
        midbrain_tmp.append(bs[1])
    else:
        midbrain_tmp.append(bs[0])

    if len(lcc) != 1:
        midbrain_tmp.append(lcc[1])
    else:
        midbrain_tmp.append(lcc[0])

    if len(lcw) != 1:
        midbrain_tmp.append(lcw[1])
    else:
        midbrain_tmp.append(lcw[0])

    if len(rcc) != 1:
        midbrain_tmp.append(rcc[1])
    else:
        midbrain_tmp.append(rcc[0])

    if len(rcw) != 1:
        midbrain_tmp.append(rcw[1])
    else:
        midbrain_tmp.append(rcw[0])

    for x in range(len(midbrain_tmp)):
        for y in range(len(midbrain_tmp[x])):
            midbrain.append(midbrain_tmp[x][y])

    ## Generate gray matter coordinates into a variable

    gm = []
    lhcoord = get_label_coord(gm_lh_labels,n_gm_lh_labels)
    rhcoord = get_label_coord(gm_rh_labels,n_gm_rh_labels)
    lhcoord_seg = get_seg_coord(get_rs_coord(lhcoord,gm_affine))
    rhcoord_seg = get_seg_coord(get_rs_coord(rhcoord,gm_affine))

    for x in range(len(lhcoord_seg[0])):
        gm.append(lhcoord_seg[0][x])
    for x in range(len(rhcoord_seg[0])):
        gm.append(rhcoord_seg[0][x])

    ## Generate ventricular coordinates into a variable

    ventricles_tmp = []
    ventricles = []
    vlhcoord = get_label_coord(seg_llv_labels,n_seg_llv_labels)
    vrhcoord = get_label_coord(seg_rlv_labels,n_seg_rlv_labels)

    if len(vlhcoord) != 1:
        ventricles_tmp.append(vlhcoord[1])
    else:
        ventricles_tmp.append(vlhcoord[0])
    if len(vrhcoord) != 1:
        ventricles_tmp.append(vrhcoord[1])
    else:
        ventricles_tmp.append(vrhcoord[0])

    for x in range(len(ventricles_tmp)):
        for y in range(len(ventricles_tmp[x])):
            ventricles.append(ventricles_tmp[x][y])

    ## Generate lesion coordinates into a variable

    lesions_les = get_label_coord(les_labels,n_les_labels)
    lesions = get_seg_coord(get_rs_coord(lesions_les,les_affine))
    lesions_seg = lesions

    ## Find centers of lesion coordinates

    les_averages = []
    for x in range(len(lesions)):
        les_averages.append(average_func(lesions[x]))

    ## Classifying lesions by lowest Euclidean distance

    min_val_jux = 2.2    #threshold for juxtacortical lesion's distance to gray matter
    min_val_per = 7.5    #threshold for periventricular lesion's distance to ventricle
    min_val_mb = 6.0     #threshold for infratentorial lesion's distance to midbrain
    les_type = []
    dis_mb = []
    dis_gm = []
    dis_v = []

    for x in range(len(les_averages)):
        gm_min = np.min(dist_det(gm,les_averages[x]))
        mb_min = np.min(dist_det(midbrain,les_averages[x]))
        v_min = np.min(dist_det(ventricles,les_averages[x]))
        sec_v_min = np.min(dist_det(lesions[x],average_func(ventricles)))
        if gm_min <= min_val_jux and gm_min <= mb_min and gm_min <= v_min:
            lesion_type = "juxtacortical"
        elif mb_min <= min_val_mb:
            lesion_type = "infratentorial"
        elif v_min <= min_val_per:
            lesion_type = "periventricular"
        elif v_min >= 13 and gm_min >= 6.5 and len(lesions[x]) >= 75 and sec_v_min <= 60:
            lesion_type = "periventricular"
        else:
            lesion_type = "subcortical"
        print "Lesion", x+1, "analysis complete"
        les_type.append(lesion_type)
        dis_gm.append(round(gm_min,decval))
        dis_mb.append(round(mb_min,decval))
        dis_v.append(round(v_min,decval))


    ## Append results to running list 

    results = []
    sub_count, inf_count, jux_count, per_count, err_count = 0, 0, 0, 0, 0
    for count in range(len(les_type)):
        results.append([count+1, les_type[count]])
        if les_type[count] == "subcortical":
            sub_count += 1
        elif les_type[count] == "infratentorial":
            inf_count += 1
        elif les_type[count] == "juxtacortical":
            jux_count += 1
        elif les_type[count] == "periventricular":
            per_count += 1
    for count in range(len(lesions)):
        all_sub_results["mseID"].append(subjects[numsubjects])
        all_sub_results["total number of lesions"].append(len(lesions))
        all_sub_results["subcortical lesions"].append(sub_count)
        all_sub_results["juxtacortical lesions"].append(jux_count)
        all_sub_results["periventricular lesions"].append(per_count)
        all_sub_results["infratentorial lesions"].append(inf_count)
        all_sub_results["lesion"].append(count+1)
        all_sub_results["type"].append(les_type[count])
        all_sub_results["center coordinates"].append(str(les_averages[count]))
        all_sub_results["volume"].append(len(lesions[count]))
        all_sub_results["distance from midbrain"].append(dis_mb[count])
        all_sub_results["distance from ventricles"].append(dis_v[count])
        all_sub_results["distance from gray matter"].append(dis_gm[count])
    print subjects[numsubjects], "data input complete."; print

## convert to csv file

sub_results = pd.DataFrame(all_sub_results,columns=["mseID", 
                                                    "total number of lesions", 
                                                    "subcortical lesions", 
                                                    "juxtacortical lesions", 
                                                    "periventricular lesions", 
                                                    "infratentorial lesions", 
                                                    "lesion", 
                                                    "type", 
                                                    "center coordinates", 
                                                    "volume",
                                                    "distance from midbrain",
                                                    "distance from ventricles",
                                                    "distance from gray matter"])
sub_results.to_csv('/data/henry1/mahamber/tmp_lesion_info.csv')
print sub_results; print
for x in range(len(problems)):
    print[problems[x]]


# In[ ]:




# In[ ]:



