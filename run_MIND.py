import sys, os
from os.path import exists
import numpy as np 
import nibabel as nib 
import scipy.io
sys.path.insert(1, '/home/riccardo/MIND/') # controlla che prende la branc MIND_dev
from MIND import compute_MIND
import glob 
import xarray as xr
import numpy as np
import pandas as pd
import pickle
from itertools import combinations
import pickle
from nipype.interfaces.freesurfer import SampleToSurface, ApplyVolTransform  # this function is a wrapper to mri_vol2surf
from utils.utils import load_img, save_img, reorient_img, parse_lut, str_to_bool, LoguruStreamHandler # type: ignore

# when you open a new terminal remember to export SUBJECTS_DIR 

def get_group(subj_id: str) -> str:
    if subj_id.startswith("sub-HC"):
        return "HC"
    elif subj_id.startswith("sub-PSP"):
        return "PSP"
    elif subj_id.startswith("sub-PD"):
        return "PD"
    elif subj_id.startswith("sub-"):  # plain numbers = MSA
        return "MSA"
    elif subj_id[0].isdigit():
        return "HCyoung"
    else:
        raise ValueError(f"Unknown group for subject {subj_id}")
    
freesurfer_path=os.environ.get("SUBJECTS_DIR")   
base_dir ="/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer"
base_dir_young = "/mnt/NAS_RAW/HCP/structuralZipped/freesurfer8"
output_path = "/home/riccardo/codici_progetti/Salerno/structural_similarity"#"/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/structural_similarity"

session=sys.argv[1] # this is the session parameter
parcellation = sys.argv[2]  # 'aparc' or 'aparc.a2009s' fo DK and Destriuex default from FreeSurfer
TH_outlier = float(sys.argv[3]) # if you want to include all datapoint (no exlusion) use something big like 99999 greater than 9000

print(TH_outlier >= 9000)
if TH_outlier >= 9000:
    outlier_name = "outliers_not_removed"
else:
    outlier_name =f"outliers_removed"

print(outlier_name)

# match only ses-T0 and ses-T1, avoid "longses"
subjects_paths = glob.glob(f"{base_dir}/sub-*_{session}") # 1]
subjects_young = glob.glob(f"{base_dir_young}/*")
subjects_young.sort()
subjects_paths.sort() 

# remove sub-PSP026, PSP029, PSP073 for which freesurfer is not succesfull
subjects_paths = [s for s in subjects_paths if "PSP027" not in s and "PSP029" not in s and "PSP073" not in s]

#subjects_paths += subjects_young  # do not run on HCP subejcts 

# rimuovi al momento questi sani dell'HCP
subjects_paths = [s for s in subjects_paths if "308129" not in s and "449753" not in s and "561444" not in s and "573451" not in s and "590047" not in s and "604537" not in s and "618952" not in s and "663755" not in s and "692964" not in s and "698168" not in s and "749058" not in s and "753150" not in s]

#index = [idx for idx, s in enumerate(subjects_paths) if '308129' in s][0]
#print(index)

all_minds = {}
all_before = {}
all_after = {} 
all_orig = {} 

groups = {"HC": [], "PD": [], "PSP": [], "MSA": [], "HCyoung": []}
features = ['CT','MC','Vol','SD','SA'] 
feature_name = ['TBM']
feature_subsets =  [features[:i] + features[i+1:] for i in range(len(features))]

# Add the full feature set as the last case
feature_subsets.append(features)
feature_subsets.reverse()  # to start from all features at first 
#subjects_paths.reverse()

for idx, subset in enumerate(feature_subsets):

    if idx > 0:
        continue  # skip all but the first iteration. remove this part of the cod eif you want to perform the leave on out 

    print(subset)
    selected_elements = [subjects_paths[i] for i in [15,201] if i < len(subjects_paths)] # 45, 88, 89, 124, 125, 199, 200, 279, 280
    #subjects_paths.reverse()
    #print(subjects_paths)

    for subj_path in subjects_paths:
        print(subj_path)

        subj_id = os.path.basename(subj_path).split("_")[0]  # e.g. "sub-01" or "sub-HC01"
        if str(subj_id).isdigit():
            subj_id = str(subj_id)

        group = get_group(subj_id)

        print(subj_path)
        MIND , before, after , vertex_data_orig = compute_MIND(subj_path, subset, parcellation, filter_vertices = True, resample = False, TH_outlier = TH_outlier) # 5th argument is resample and must be TRUE when doing with univariate approach see git hub page of MIND, default is False
        #print(type(before), type(MIND), before.shape)
        # Store subject-level
        subj_str = str(subj_id)
        if subj_str[0].isdigit():
            print("########### HC young ###########")
            safe_key = "sub-" + subj_str
        else: 
            safe_key = subj_id

        all_minds[safe_key] = {'SubjID': subj_str , 'matrix': MIND, 'group': group}
        all_before[safe_key] = {'SubjID': subj_str , 'matrix': before, 'group': group}
        all_after[safe_key] = {'SubjID': subj_str , 'matrix': after, 'group': group}
        all_orig[safe_key] = {'SubjID': subj_str , 'matrix': vertex_data_orig, 'group': group}

        print(f"end of subject {subj_id}")


        # Save both dictionaries together in one file
    with open(f"{output_path}/MIND_networks_results/MIND_networks_ses-T0_{'_'.join(subset)}_{parcellation}_{outlier_name}_{TH_outlier}.pkl", "wb") as f:
        pickle.dump({"all_minds": all_minds}, f)
    with open(f"{output_path}/distributions/before_outlier_removal_ses-T0_{'_'.join(subset)}_{parcellation}_{outlier_name}_{TH_outlier}.pkl", "wb") as f:
        pickle.dump({"all_before": all_before}, f)
    with open(f"{output_path}/distributions/after_outlier_removal_ses-T0_{'_'.join(subset)}_{parcellation}_{outlier_name}_{TH_outlier}.pkl", "wb") as f:
        pickle.dump({"all_after": all_after}, f)
    with open(f"{output_path}/distributions/vertex_data_orig_ses-T0_{'_'.join(subset)}_{parcellation}_{outlier_name}_{TH_outlier}.pkl", "wb") as f:
        pickle.dump({"all_orig": all_orig}, f)
        
 