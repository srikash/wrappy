#!/bin/env python
"""
Python scripts for commonly used pipelines
"""
import os
import sys
from . import misc, fsl
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp

def epi_proc(input_file_path, ref_file_path=None,tr=1.0,topup_prefix=None,topup_parameters=None,topup_in_index=1.0):
    """Motion + distortion correction workflow

    Args:
        input_file_path (PosixPath): Path to input file.
        ref_file_path (PosixPath, optional): Path to reference file. Defaults to None.
        tr (float, optional): Volume TR. Defaults to 1.0
        topup_prefix (str, optional): Path to topup output prefix. Defaults to None.
        topup_parameters (PosixPath, optional): Path to topup parameters file. Defaults to None.
        topup_in_index (int, optional): Input order. Defaults to 1.0.
    """
    out_prefix = misc.get_prefix(input_file_path)
    
    fsl.fsl_split(input_file_path=input_file_path)

    in_splitlist = list(out_prefix.glob("vol_*[!fslWarped].nii.gz"))
    in_splitlist.sort()  # To make sure volumes are ordered correctly
    
    if ref_file_path is None:
        misc.make_identity_matrix(out=out_prefix.as_posix() + "/vol_0000_fsl.mat")
        fsl.fsl_realign(input_file_list=in_splitlist)
        out_splitlist_coreg = list(out_prefix.glob("vol_*_fslWarped.nii.gz"))
        out_splitlist_coreg.sort()  # To make sure volumes are ordered correctly
        out_mergelist_0 = misc.list_to_string(in_list=out_splitlist_coreg)
        out_mergelist = in_splitlist[0].as_posix() + " " + out_mergelist_0
    else:
        fsl.fsl_realign(input_file_list=in_splitlist,fixed_file_path=ref_file_path)
        out_splitlist_coreg = list(out_prefix.glob("vol_*_fslWarped.nii.gz"))
        out_splitlist_coreg.sort()  # To make sure volumes are ordered correctly
        out_mergelist = misc.list_to_string(in_list=out_splitlist_coreg)
        
    fsl.fsl_merge(input_file_list=out_mergelist,out_file_path=Path(out_prefix.as_posix()+"_fslMoCorr.nii.gz"),tr=tr)

    misc.remove_files(src=out_prefix,pattern="vol_*.nii.gz")
    
    out_tmean = fsl.fsl_Tmean(input_file_path=Path(out_prefix.as_posix() + "_fslMoCorr.nii.gz"))

    if topup_prefix and topup_parameters is not None:
        fsl.fsl_applytopup(input_file_path=Path(out_prefix.as_posix()+"_fslMoCorr.nii.gz"),
                           topup_parameters=topup_parameters,
                           topup_output_prefix=topup_prefix,
                           in_index=topup_in_index)
    
    return out_tmean

def prep_topup(pe_file_path, ope_file_path,datain_file, config_file=None):
    """ Prepare topup
    
    Args:
        pe_file_path (PosixPath): Path to PE file.
        ope_file_path (PosixPath): Path to oPE file.
        topup_parameters (PosixPath): Path to topup parameters file
        topup_config (PosixPath, optional): Path to topup configuration. Defaults to topup default.
    """
    out_prefix = misc.get_prefix(pe_file_path)
    out_string = str(out_prefix.as_posix()).replace("_dir-PE_run-01_fslMoCorr_fslTmean","")
    
    moving_file_out = fsl.fsl_register(fixed_file_path=pe_file_path, moving_file_path=ope_file_path)
    merge_list = [pe_file_path]
    merge_list.append(Path(moving_file_out))
    
    fsl.fsl_merge(input_file_list=misc.list_to_string(in_list=merge_list),out_file_path=Path(out_string + "_infile.nii.gz"))

    if config_file is not None:
        fsl.fsl_topup(input_file_path=Path(out_string + "_infile.nii.gz"), topup_parameters=datain_file, topup_config=config_file)    
    else:
        fsl.fsl_topup(input_file_path=Path(out_string + "_infile.nii.gz"), topup_parameters=datain_file)

def run_oxasl(asl_file_path,mzero_file_path,mask_file_path,asl_data="diff"):
    """[summary]

    Args:
        asl_file_path (PosixPath): Path to ASL data file
        mzero_file_path (PosixPath): Path to M0 data file
        mask_file_path (PosixPath): Path to M0 brainmask
        asl_data (str, optional): "tc" or "diff". Defaults to "diff".
    """
    
    out_prefix = misc.get_prefix(asl_file_path)
    out_folder = Path(out_prefix.as_posix() + "_oxasl")
    
    main_cmd = "oxasl --overwrite --te=14 --mode=longtr --tr=20.0" + " " + \
        "--calib-aslreg --calib-alpha=0.95 --calib-method=voxelwise" + " " + \
        "--t1b=2.2 --t1=1.8 --fixbolus --spatial-off --slicedt=0.03197" + " " + \
        "--iaf=" + asl_data + " " + "--ibf=rpt --ntis=1 --bolus=0.7 --tis=1.8" + " " + \
        "--calib=" + mzero_file_path.as_posix() + " " + \
        "--mask=" + mask_file_path.as_posix() + " " + \
        "--asldata=" + asl_file_path.as_posix() + " " + \
        "--output=" + out_folder.as_posix()
    misc.exec_shell(cmd=main_cmd)

    shutil.copyfile(out_folder.as_posix()+"/output/native/perfusion.nii.gz",out_prefix.as_posix()+"_perfusion.nii.gz")
    shutil.copyfile(out_folder.as_posix()+"/output/native/perfusion_calib.nii.gz",out_prefix.as_posix()+"_perfusion_calib.nii.gz")