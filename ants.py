#!/bin/env python
"""
Python scripts for ANTs command-line utilities
"""
import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp

def resample_image(input_file_path, resolution=0.3, verbose=False, sim=False):
    """Resamples an image to a given resolution isotropic.

    Args:
        input_file_path (PosixPath): Path to input file
        resolution (float, optional): Resolution to resample to. Defaults to 0.3.
        verbose (bool, optional): Prints loads of stuff to terminal. Defaults to False.
    """
    def get_prefix(src):
        """Removes prefix from Path file

        Args:
            src (PosixPath): path to the input file

        Returns:
            PosixPath: path to the input file without extension
        """
        check_suffix = src.suffixes
        if len(check_suffix) == 1:
            src_prefix = src.with_suffix('')
        elif len(check_suffix) == 2:
            src_prefix = src.with_suffix('').with_suffix('')
        return src_prefix

    if type(input_file_path) is str:
        print("Requires Input as PosixPath")

    out_file_prefix = get_prefix(input_file_path)
    out_file_path = Path(out_file_prefix.as_posix() + "_0p3_iso.nii.gz")

    main_cmd = "${ANTSPATH}/ResampleImage " + \
        "3 " + \
        input_file_path.as_posix() + \
        " " + \
        out_file_path.as_posix() + \
        " " + \
        str(resolution)+"x"+str(resolution)+"x"+str(resolution) + \
        " " + \
        "0 4[5] 2 "

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd)
    if sim is False:
        if verbose is False:
            subprocess.run(main_cmd, shell=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elif verbose is True:
            subprocess.run(main_cmd, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")
    return out_file_path


def build_template(input_file_list, outprefix, verbose=False, sim=False):
    """Runs ANTs template construction command.

    Args:
        input_file_list (list): List of files to be made template out of
        out_prefix (str): Output name prefix
        verbose (bool, optional): Prints loads of stuff to terminal. Defaults to False.
    """
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

    if type(input_file_list) is str:
        print("Requires Input as List")

    input_file_string = list_to_string(in_list=input_file_list)

    main_cmd = "${ANTSPATH}/antsMultivariateTemplateConstruction2.sh " + \
        "-d 3 -a 0 -i 4 -k 1 -f 4x2x1 -s 2x1x0vox " + \
        "-q 30x20x4 -t SyN -m CC[2] -c 2 -j 34 " + \
        "-o " + outprefix + \
        " " + \
        input_file_string

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd)
    if sim is False:
        if verbose is False:
            subprocess.run(main_cmd, shell=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elif verbose is True:
            subprocess.run(main_cmd, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")


def apply_transform(ref_file_path, trg_file_path, mat_file_path, warp_file_path=None, verbose=False, sim=False):
    """Apply estimated transformations using ANTs

    Args:
        ref_file_path (PosixPath): Path to the reference image space
        trg_file_path (PosixPath): Path to data being transformed
        mat_file_path (PosixPath): Path to affine matrix 
        warp_file_path (PosixPath, optional): Path to warp. Defaults to None.
        verbose (bool, optional): [description]. Defaults to False.
        sim (bool, optional): [description]. Defaults to False.
    """
    def get_prefix(src):
        """Removes prefix from Path file

        Args:
            src (PosixPath): path to the input file

        Returns:
            PosixPath: path to the input file without extension
        """
        check_suffix = src.suffixes
        if len(check_suffix) == 1:
            src_prefix = src.with_suffix('')
        elif len(check_suffix) == 2:
            src_prefix = src.with_suffix('').with_suffix('')
        return src_prefix

    if type(ref_file_path) is str:
        print("Requires Input as PosixPath")
    elif type(trg_file_path) is str:
        print("Requires Input as PosixPath")

    out_file_prefix = get_prefix(trg_file_path)
    out_file_path = Path(out_file_prefix.as_posix() + "_Warped.nii.gz")

    if warp_file_path is None:
        main_cmd = "${ANTSPATH}/antsApplyTransforms " + \
            "-d 3 -u float -n BSpline[5] " + \
            "-r " + ref_file_path.as_posix() + \
            " " + \
            "-i " + trg_file_path.as_posix() + \
            " " + \
            "-o " + out_file_path.as_posix() + \
            " " + \
            "-t " + mat_file_path.as_posix()
    else:
        main_cmd = "${ANTSPATH}/antsApplyTransforms " + \
            "-d 3 -u float -n BSpline[5] " + \
            "-r " + ref_file_path.as_posix() + \
            " " + \
            "-i " + trg_file_path.as_posix() + \
            " " + \
            "-o " + out_file_path.as_posix() + \
            " " + \
            "-t " + warp_file_path.as_posix() + \
            " " + \
            "-t " + mat_file_path.as_posix()

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd)
    if sim is False:
        if verbose is False:
            subprocess.run(main_cmd, shell=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elif verbose is True:
            subprocess.run(main_cmd, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")
    return out_file_path


def average_images(input_file_list, verbose=False):
    """Averages a list of given images.

    Args:
        input_file_list (list): List of files to be averaged
        verbose (bool, optional): Prints loads of stuff to terminal. Defaults to False.
    """
    def get_prefix(src):
        """Removes prefix from Path file

        Args:
            src (PosixPath): path to the input file

        Returns:
            PosixPath: path to the input file without extension
        """
        check_suffix = src.suffixes
        if len(check_suffix) == 1:
            src_prefix = src.with_suffix('')
        elif len(check_suffix) == 2:
            src_prefix = src.with_suffix('').with_suffix('')
        return src_prefix

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

    if type(input_file_list) is str:
        print("Requires Input as List")

    input_file_string = list_to_string(in_list=input_file_list)

    out_file_prefix = get_prefix(input_file_list[0])
    out_file_path = Path(out_file_prefix.as_posix() + "_Avg.nii.gz")

    main_cmd = "${ANTSPATH}/AverageImages " + \
        "3 " + \
        out_file_path.as_posix() + \
        " " + \
        "0" + \
        " " + \
        input_file_string

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd)
    if verbose is False:
        subprocess.run(main_cmd, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif verbose is True:
        subprocess.run(main_cmd, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")
    return out_file_path