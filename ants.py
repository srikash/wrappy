#!/bin/env python
"""
Python scripts for ANTs command-line utilities
"""
import os
import sys
import misc
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp

def ants_resample(input_file_path, resolution=0.3, interp="spline",verbose=False, sim=False):
    """Resamples an image to a given isotropic resolution.

    Args:
        input_file_path (PosixPath): Path to input file
        resolution (float, optional): Resolution to resample to. Defaults to 0.3.
        interp (str, optional): "spline" or "sinc". Defaults to "spline". 
        verbose (bool, optional): Prints loads of stuff to terminal. Defaults to False.
    """
    out_file_prefix = misc.get_prefix(input_file_path)
    res_str = str(resolution).replace(".","p")
    out_file_path = Path(str(out_file_prefix.as_posix() + "_"+res_str+"_"+"iso.nii.gz")) 
    
    if interp == "spline":
        interpolator = str("4[5]")
    elif interp == "sinc":
        interpolator = str("3['l']")

    main_cmd = "ResampleImage 3 " + \
        input_file_path.as_posix() + " " + \
        out_file_path.as_posix() + " " + \
        str(resolution)+"x"+str(resolution)+"x"+str(resolution) + " " + \
        "0 " + interpolator + " 2"

    misc.exec_shell(cmd=main_cmd)
    return out_file_path


def ants_build_template(input_file_list, out_prefix, iters=4,verbose=False, sim=False):
    """Runs ANTs template construction command.

    Args:
        input_file_list (list): List of files to be made template out of
        out_prefix (str): Output name prefix
        iters (int, optional): Number of iterations. Defaults to 4.
        verbose (bool, optional): Print verbose output. Defaults to False.
    """
    num_cores = mp.cpu_count()
        
    if type(input_file_list) is str:
        print("Requires Input as List")

    input_file_string = misc.list_to_string(in_list=input_file_list)

    main_cmd = "antsMultivariateTemplateConstruction2.sh " + \
        "-d 3 -a 0 " + \
        "-k 1 -f 4x2x1 -s 4x2x0vox " + \
        "-q 50x20x5 -t SyN -m CC[2] -c 2 " + \
        "-j " + num_cores + " " + \
        "-i " + iters + " " + \
        "-o " + out_prefix + " " + \
        input_file_string

    misc.exec_shell(cmd=main_cmd)

def ants_reslice(fixed_file_path, moving_file_path, mat_file_path, warp_file_path=None, interp="spline", verbose=False):
    """Apply transformations using antsApplyTransforms

    Args:
        fixed_file_path (PosixPath): Path to the fixed image
        moving_file_path (PosixPath): Path to the moving image
        mat_file_path (PosixPath): Path to the affine matrix
        warp_file_path (PosixPath, optional): Path to the warp. Defaults to None.
        interp (str, optional): "spline" or "sinc". Defaults to "spline".
        verbose (bool, optional): Print verbose output. Defaults to False.
    
    TO-DO: Accomodate many mats and warps as lists 
    """
    if type(fixed_file_path) is str:
        fixed_file_path=Path(fixed_file_path)
    elif type(moving_file_path) is str:
        moving_file_path=Path(moving_file_path)

    out_file_prefix = misc.get_prefix(moving_file_path)
    out_file_path = Path(out_file_prefix.as_posix() + "_antsWarped.nii.gz")

    if interp == "spline":
        interpolator = str("BSpline[5]")
    elif interp == "sinc":
        interpolator = str("LanczosWindowedSinc")
    
    if warp_file_path is None:
        main_cmd = "antsApplyTransforms " + \
            "-d 3 -u float " + \
            "-n " + interpolator + " " + \
            "-r " + fixed_file_path.as_posix() + " " + \
            "-i " + moving_file_path.as_posix() + " " + \
            "-o " + out_file_path.as_posix() + " " + \
            "-t " + mat_file_path.as_posix()
        misc.exec_shell(cmd=main_cmd)
    else:
        main_cmd = "antsApplyTransforms " + \
            "-d 3 -u float " + \
            "-n " + interpolator + " " + \
            "-r " + fixed_file_path.as_posix() + " " + \
            "-i " + moving_file_path.as_posix() + " " + \
            "-o " + out_file_path.as_posix() + " " + \
            "-t " + warp_file_path.as_posix() + " " + \
            "-t " + mat_file_path.as_posix()
        misc.exec_shell(cmd=main_cmd)
    return out_file_path


def ants_average(input_file_list, norm=0,verbose=False):
    """Calculate average using AverageImages

    Args:
        input_file_list (list): List of files to be averaged
        norm (int, optional):   0 --> no normalisation
                                1 --> divide by mean + sharpening
                                2 --> divide by mean 
                                Defaults to 0.
        verbose (bool, optional): Print verbose output. Defaults to False.
    """

    if type(input_file_list) is str:
        print("Requires Input as List")

    input_file_string = misc.list_to_string(in_list=input_file_list)

    out_file_prefix = misc.get_prefix(input_file_list[0])
    out_file_path = Path(out_file_prefix.as_posix() + "_antsAvg.nii.gz")

    main_cmd = "AverageImages 3 " + \
        out_file_path.as_posix() + " " + \
        int(norm) + " " + \
        input_file_string

    misc.exec_shell(cmd=main_cmd)

def ants_imath(input_file_path,operation=None,parameters=None):
    """Image operations using ImageMath

    Args:
        input_file_path (PosixPath): Path to the input file.
        operation (str, optional):  "anismooth"
                                    "close"
                                    "contours"
                                    "dilate"
                                    "erode"
                                    "fillholes"
                                    "largest"
                                    "mask"
                                    "negative"
                                    "open"
                                    "projmax"
                                    "projmin"
                                    "projsum"
                                    "rescale"                                    
                                    "truncate"                                    
                                    Defaults to None.
        parameters (str, optional): "anismooth" --> iters conductance
                                                    10    0.25
                                    "close"     --> voxels
                                                    2
                                    "contours"
                                    "dilate"    --> voxels
                                                    2
                                    "erode"     --> voxels
                                                    2
                                    "fillholes"
                                    "largest"
                                    "mask"      --> % of mean
                                                    1.2
                                    "negative"
                                    "open"      --> voxels
                                                    2
                                    "projmax"   --> axis
                                                    3
                                    "projmin"   --> axis
                                                    3
                                    "projsum"   --> axis
                                                    3
                                    "rescale"   --> min max
                                                    0   4000                                    
                                    "truncate"  --> %min %max
                                                    0.01 0.99                                    
                                    Defaults to None.
    """
    if type(input_file_path) is str:
        input_file_path=Path(input_file_path)

    out_file_prefix = misc.get_prefix(input_file_path)
    
    if operation == "anismooth":
        operation = "PeronaMalik"
        if parameters is None:
            parameters="10 0.25"
    elif operation == "close":
        operation = "MC"
        if parameters is None:
            parameters="2"               
    elif operation == "contours":
        operation = "ExtractContours"
        if parameters is None:
            parameters="1"
    elif operation == "dilate":
        operation = "MD"
        if parameters is None:
            parameters="2"
    elif operation == "erode":
        operation = "ME"
        if parameters is None:
            parameters="2"
    elif operation == "fillholes":
        operation = "FillHoles"
        if parameters is None:
            parameters="2"
    elif operation == "largest":
        operation = "GetLargestComponent"
    elif operation == "mask":
        operation = "ThresholdAtMean"
        if parameters is None:
            parameters="1.2"
    elif operation == "negative":
        operation = "Neg"
    elif operation == "open":
        operation = "MO"
        if parameters is None:
            parameters="2"
    elif operation == "projmax":
        operation = "Project"
        if parameters is None:
            parameters="3"
        parameters = parameters + " 1"
    elif operation == "projmin":
        operation = "Project"
        if parameters is None:
            parameters="3"
        parameters = parameters + " 2"
    elif operation == "projsum":
        operation = "Project"
        if parameters is None:
            parameters="3"
        parameters = parameters + " 0"
    elif operation == "rescale":
        operation = "RescaleImage"
        if parameters is None:
            parameters="0 4000"
    elif operation == "truncate":
        operation = "TruncateImageIntensity"
        parameters = "0.01 0.99"

    out_file_path=Path(out_file_prefix.as_posix()+"_"+operation+".nii.gz")
    
    main_cmd = "ImageMath 3 " + \
        out_file_path.as_posix() + " " + \
        operation + " " + \
        input_file_path.as_posix() + " " + \
        parameters
    
    misc.exec_shell(cmd=main_cmd)    

    return out_file_path

def ants_n4(input_file_path,mask_file_path=None,verbose=False):
    """[summary]

    Args:
        input_file_path (PosixPath): Path to the input file.
        mask_file_path (PosixPath, optional): Path to the mask file. Defaults to None.
        verbose (bool, optional): Print verbose output. Defaults to False.
    """
    
    if type(input_file_path) is str:
        input_file_path=Path(input_file_path)

    out_file_prefix = misc.get_prefix(input_file_path)
    n4_out_file_path = Path(out_file_prefix.as_posix() + "_antsN4corr.nii.gz")
    bias_out_file_path = Path(out_file_prefix.as_posix() + "_antsN4bias.nii.gz")
    
    trunc_in_file_path = ants_imath(input_file_path,operation="truncate")
    
    if mask_file_path is not None:
        main_cmd = "N4BiasFieldCorrection " + \
                "-d 3 -s 2 " + \
                "-x " + mask_file_path.as_posix() + " " + \
                "-i " + trunc_in_file_path.as_posix() + " " + \
                "-o [" + n4_out_file_path.as_posix() + \
                    ", " + bias_out_file_path.as_posix() + "]"

        misc.exec_shell(cmd=main_cmd)
    else:
        main_cmd = "N4BiasFieldCorrection " + \
                "-d 3 -s 2 " + \
                "-i " + input_file_path.as_posix() + " " + \
                "-o [" + n4_out_file_path.as_posix() + \
                    ", " + bias_out_file_path.as_posix() + "]"
        misc.exec_shell(cmd=main_cmd)

def ants_kmeans(input_file_path,mask_file_path,classes=3,verbose=False):
    """Do K-means segmentation using ANTs

    Args:
        input_file_path (PosixPath): Path to the input file
        mask_file_path (PosixPath): Path to the mask file.
        classes (int, optional): [description]. Defaults to 3.
        verbose (bool, optional): [description]. Defaults to False.
    """
    if type(input_file_path) is str:
        input_file_path=Path(input_file_path)

    out_file_prefix = misc.get_prefix(input_file_path)
    classes_out_file_path = Path(out_file_prefix.as_posix() + "_antsClasses.nii.gz")
    posteriors_out_file_path = Path(out_file_prefix.as_posix() + "_antsPosteriors.nii.gz")
    
    main_cmd = "Atropos -d 3 -c 10" + " " + \
        "-i KMeans[" + str(classes) +"]" + " " + \
        "-a " + input_file_path.as_posix() + " " + \
        "-x " + mask_file_path.as_posix() + " " + \
        "-o [" + classes_out_file_path.as_posix() + \
            ", " + posteriors_out_file_path.as_posix() + "]"
    misc.exec_shell(cmd=main_cmd)

def ants_otsu(input_file_path,mask_file_path,classes=3,verbose=False):
    """Do Otsu segmentation using ANTs

    Args:
        input_file_path (PosixPath): Path to the input file
        mask_file_path (PosixPath): Path to the mask file.
        classes (int, optional): [description]. Defaults to 3.
        verbose (bool, optional): [description]. Defaults to False.
    """
    if type(input_file_path) is str:
        input_file_path=Path(input_file_path)

    out_file_prefix = misc.get_prefix(input_file_path)
    classes_out_file_path = Path(out_file_prefix.as_posix() + "_antsClasses.nii.gz")
    posteriors_out_file_path = Path(out_file_prefix.as_posix() + "_antsPosteriors.nii.gz")
    
    main_cmd = "Atropos -d 3 -c 10" + " " + \
        "-i Otsu[" + str(classes) +"]" + " " + \
        "-a " + input_file_path.as_posix() + " " + \
        "-x " + mask_file_path.as_posix() + " " + \
        "-o [" + classes_out_file_path.as_posix() + \
            ", " + posteriors_out_file_path.as_posix() + "]"
    misc.exec_shell(cmd=main_cmd)