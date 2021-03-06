#!/bin/env python
"""
Python scripts for commonly used pipelines
"""
import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp

def epi_proc(src, ref=None):
    """[summary]

    Args:
        src ([type]): [description]
        to_zero (bool, optional): [description]. Defaults to False.

    Returns:
        [type]: [description]
    """
    print(" ++ Making directory")
    print(" ++ Splitting 4D file")
    split4d(src)
    src_prefix = get_prefix(src)
    src_splitlist = list(src_prefix.glob("vol_*[!coreg].nii.gz"))
    src_splitlist.sort()  # To make sure volumes are ordered correctly
    if ref is None:
        make_identity_matrix(out=str(src_prefix)+"/vol_0000_coreg.mat")
        print(" ++ Running registrations in parallel")
        mcflirt(file_list=src_splitlist)
    src_splitlist_coreg = list(src_prefix.glob("vol_*_coreg.nii.gz"))
    src_splitlist_coreg.sort()  # To make sure volumes are ordered correctly
    src_mergelist = list_to_string(in_list=src_splitlist_coreg)
    src_mergelist = src_splitlist[0].as_posix() + " " + src_mergelist
    print(" ++ Merging to 4D file")
    merge3d(src=src_mergelist,
            out=str(src_prefix)+"_moCorr.nii.gz",
            tr=20.0)
    print(" ++ Cleaning up")
    remove_files(src=src_prefix,
                 pattern="vol_*.nii.gz")
    print(" ++ Calculating Tmean of realigned volumes")
    src_tmean = calc_Tmean(src=str(src_prefix)+"_moCorr.nii.gz")
    return src_tmean
