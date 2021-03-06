#!/bin/env python
"""
Python scripts for miscellaneous utility
"""
import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path
from scipy import ndimage
import multiprocessing as mp
import matplotlib.pyplot as plt
from skimage import transform, util

def exec_parashell(cmd_list):
    num_cores = mp.cpu_count()
    processing_pool = mp.Pool(num_cores)
    processing_pool.map(os.system, (cmd for cmd in cmd_list))

def get_prefix(src, verbose=False):
    """Removes extension and gets the prefix from Path

    Args:
        src (PosixPath): path to the input file

    Returns:
        PosixPath: path to the input file without extension
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Checking File Extension")
        print("%---------------------------------------%")
    check_suffix = src.suffixes
    if len(check_suffix) == 1:
        if verbose is True:
            print("++++ File has a NIFTI suffix.")
            print("%---------------------------------------%")
        src_prefix = src.with_suffix('')
    elif len(check_suffix) == 2:
        if verbose is True:
            print("++++ File has a NIFTI_GZ suffix.")
            print("%---------------------------------------%")
        src_prefix = src.with_suffix('').with_suffix('')
    return src_prefix

def make_dir(src, verbose=False):
    """Makes a directory with the name of the file

    Args:
        src (PosixPath): path to the input file
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Making Directory")
        print("%---------------------------------------%")
    if len(src.suffixes) != 0:
        newfolder_path = get_prefix(src)
    else:
        newfolder_path = src
    if verbose is True:
        print("++++ Directory Name : ")
        print("     "+str(newfolder_path))
        print("%---------------------------------------%")
    newfolder_path.mkdir(parents=True, exist_ok=True)


def remove_dir(src, verbose=False):
    """Remove a directory with the name of the file

    Args:
        src (PosixPath): path to the input file
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Removing Directory")
        print("%---------------------------------------%")
    if len(src.suffixes) != 0:
        newfolder_path = get_prefix(src)
    else:
        newfolder_path = src
    if verbose is True:
        print("++++ Directory Name : ")
        print("     "+str(newfolder_path))
        print("%---------------------------------------%")
    shutil.rmtree(newfolder_path)


def list_to_string(in_list):
    """Convert list to string, by joining all item in list.

    Args:
        in_list (list): input list
        seperator (str, optional): delimiter. Defaults to ' '.

    Returns:
        string: Returns the concatenated string
    """
    full_str = ' '.join([str(elem.as_posix()) for elem in in_list])
    return full_str


def remove_files(src, pattern):
    """Removes files matching a certain pattern from src

    Args:
        src (PosixPath): path to directory
        pattern (string): pattern in files
    """
    for file in src.glob(pattern):
        file.unlink()


def make_identity_matrix(out="identity.mat"):
    """Creates an identity matrix for purposes
    Args:
        out (str, optional): path to output filename. Defaults to "identity.mat".
    """
    # Fix the Posix path in case
    if type(out) is not str:
        out = out.as_posix()
    np.savetxt(out, np.identity(4, dtype=int), fmt="%i", delimiter=" ")

def make_checkerboard(size=16):
    """[summary]

    Args:
        size (int, optional): size of checkered square. Defaults to 16.

    Returns:
        [type]: [description]
    """

    checkerboard = np.zeros((size,size),dtype=int)
    checkerboard[1::2,::2] = 1
    checkerboard[::2,1::2] = 1
    checker_size = checkerboard.shape[0]
    
    return checkerboard,checker_size

def make_montage(input_file_path,slice_range=np.arange(13,22)):
    input_nii=nib.load(input_file_path.as_posix())
    input_nii_data = input_nii.get_fdata()
    
    data_subset=np.zeros([len(slice_range),input_nii_data.shape[0],input_nii_data.shape[1]])
    
    for idx in range(0,len(slice_range)):
        data_subset[idx,:,:]=ndimage.rotate(input_nii_data[:,:,slice_range[idx]],90)
    
    data_subset_montage = util.montage(data_subset)
    data_subset_size = data_subset.shape
    return data_subset_montage, data_subset_size

def plot_montage(input_image,colour,title=):
    fig, ax = plt.subplots(figsize=(50,15))
    ax.imshow(input_image,cmap=colour)