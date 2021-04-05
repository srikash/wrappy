#!/bin/env python
"""
Python scripts for AFNI command-line utilities
"""
import os
import sys
import shutil
from . import misc
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp


def afni_getvols(input_file_path, list_of_volumes):
    """Use AFNI to extract volumes from a timeseries

    Args:
        input_file_path (PosixPath): Path to input file
        list_of_volumes (list, int): list of timeries sub-volumes (Note: zero indexing)
    """

    if list_of_volumes is list:
        string_of_subvols = misc.list_to_string(in_list=list_of_volumes)
    else:
        string_of_subvols = list_of_volumes

    out_file_prefix = misc.get_prefix(input_file_path)

    main_cmd = "3dcalc" + " " + \
        "-overwrite" + " " + \
        "-prefix " + out_file_prefix.as_posix() + "_subvols.nii.gz" + " " + \
        "-expr a" + " " + \
        "-a " + input_file_path.as_posix() + "[" + string_of_subvols + "]"

    misc.exec_shell(cmd=main_cmd)


def afni_Tstats(input_file_path, operation="mean"):
    """ Temporal statistics using 3dTstat

    Args:
        input_file_path (PosixPath): Input file path
        operation (str, optional): Calculation to perform. Defaults to "mean".
                                   (e.g. sum, mean, stdev, tsnr, cvarinv, 
                                    skewness, kurtosis, autocorr, tdiff etc.)

    Returns:
        PosixPath: Output file path
    """

    out_file_prefix = misc.get_prefix(input_file_path)
    out_string = out_file_prefix.as_posix() + "_Tmean.nii.gz"

    main_cmd = "3dTstat" + " " + \
        "-overwrite" + " " + \
        "-prefix " + out_string + " " + \
        "-" + operation + " " + input_file_path.as_posix()

    misc.exec_shell(cmd=main_cmd)
    return Path(out_string)


def afni_volreg(input_file_path):
    """ Run 3dvolreg with some default parameters 

    Args:
        input_file_path (PosixPath): Input file path
    """

    out_file_prefix = misc.get_prefix(input_file_path)

    main_cmd = "3dvolreg" + " " + \
        "-overwrite -verbose -nomaxdisp" + " " + \
        "-nocoarserot -twopass -final heptic" + " " + \
        "-maxite 42 -noclip -zpad 2 -wtinp " + " " + \
        "-prefix " + out_file_prefix.as_posix()+"_volReg.nii.gz" + " " + \
        "-1Dfile " + out_file_prefix.as_posix()+"_volReg.params" + " " + \
        input_file_path.as_posix()

    misc.exec_shell(cmd=main_cmd)
