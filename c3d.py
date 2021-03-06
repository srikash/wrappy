#!/bin/env python
"""
Python scripts for C3D and greedy
"""
import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp

def c3d_convert(ref, src, mat):
    """Convert FSL matrix to ITK

    Args:
        ref (PosixPath): path to reference file
        src (PosixPath): path to source file
        mat (PosixPath): path to mat to be inverted
        invert (float, optional): invert matrix. Defaults to 0.0.
    """
    if type(src) is str:
        src = Path(src)
    out = get_prefix(src)

    main_cmd = "c3d_affine_tool " + \
        "-ref " + \
        str(ref) + \
        " " + \
        "-src " + \
        str(src) + \
        " " + \
        str(mat) + \
        " " + \
        "-fsl2ras " + \
        "-oitk " + \
        str(out)+"_itk.txt"
    print(main_cmd)
    subprocess.run(main_cmd, shell=True)

    inv_cmd = "c3d_affine_tool " + \
        "-ref " + \
        str(ref) + \
        " " + \
        "-src " + \
        str(src) + \
        " " + \
        str(mat) + \
        " " + \
        "-fsl2ras " + \
        "-inv " + \
        "-oitk " + \
        str(out)+"_itk_inverse.txt"
    print(inv_cmd)
    subprocess.run(inv_cmd, shell=True)
    
def greedy_register(rage_file_path,tse_file_path,verbose=False):
    """Register the TSE to hires MP2RAGE using greedy

    Args:
        rage_file_path (PosixPath): Path to hires MP2RAGE file
        tse_file_path (PosixPath): Path to T2w TSE file
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

    out_file_prefix = get_prefix(tse_file_path)
    
    main_cmd_0 = "greedy -d 3 -dof 6 -a -ia-identity " + \
        "-threads 34 -m NMI -n 50x25 " + \
        "-i " + rage_file_path.as_posix() + " " + tse_file_path.as_posix() + " " + \
        "-o " + out_file_prefix.as_posix() + "_greedy.mat"

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd_0)
    if verbose is False:
        subprocess.run(main_cmd_0, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif verbose is True:
        subprocess.run(main_cmd_0, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")
    
    main_cmd_1 = "c3d_affine_tool " + \
    "-ref " + rage_file_path.as_posix() + \
    " " + \
    "-src " + tse_file_path.as_posix() + \
    " " + \
    out_file_prefix.as_posix() +"_greedy.mat" + \
    " " + \
    "-ras2fsl " + \
    "-o " + out_file_prefix.as_posix()  + "_fsl.mat"

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd_1)
    if verbose is False:
        subprocess.run(main_cmd_1, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif verbose is True:
        subprocess.run(main_cmd_1, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")
    
    main_cmd_2 = "convert_xfm " + \
    "-omat " + out_file_prefix.as_posix()  + "_fsl_inverse.mat" + \
    " " + \
    "-inverse " + out_file_prefix.as_posix() +"_fsl.mat"

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd_2)
    if verbose is False:
        subprocess.run(main_cmd_2, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif verbose is True:
        subprocess.run(main_cmd_2, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")

    main_cmd_3 = "c3d_affine_tool " + \
    "-ref " + rage_file_path.as_posix() + \
    " " + \
    "-src " + tse_file_path.as_posix() + \
    " " + \
    out_file_prefix.as_posix() +"_greedy.mat" + \
    " " + \
    "-oitk " + out_file_prefix.as_posix() + "_itk.mat"

    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Executing process                          ")
    print("% ---------------------------------------------------------------- %")
    print(" ")
    print(" ")
    print(main_cmd_3)
    if verbose is False:
        subprocess.run(main_cmd_3, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif verbose is True:
        subprocess.run(main_cmd_3, shell=True)
    print(" ")
    print(" ")
    print("% ---------------------------------------------------------------- %")
    print("                         Process Completed                          ")
    print("% ---------------------------------------------------------------- %")


def greedy_transform(ref_file_path, trg_file_path, mat_file_path, verbose=False, sim=False):
    """Apply estimated transformations using C3D

    Args:
        ref_file_path (PosixPath): Path to the reference image space
        trg_file_path (PosixPath): Path to data being transformed
        mat_file_path (PosixPath): Path to affine matrix 
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
    out_file_path = Path(out_file_prefix.as_posix() + "_c3dWarped.nii.gz")

    main_cmd = "c3d " + \
        ref_file_path.as_posix() + \
        " " + \
        trg_file_path.as_posix() + \
        " " + \
        "-interpolation Sinc " + \
        "-reslice-matrix " + mat_file_path.as_posix() + \
        " " + \
        "-o " + out_file_path.as_posix()

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
