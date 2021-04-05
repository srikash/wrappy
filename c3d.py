#!/bin/env python
"""
Python scripts for c3d and greedy
"""
import os
import sys
# sys.path.append(os.getcwd())
from . import misc
from pathlib import Path
import multiprocessing as mp


def c3d_convert(fixed_file_path, moving_file_path, mat_file_path, ras2fsl=True, ras2itk=False, fsl2ras=False, fsl2itk=False):
    """Wrapper for c3d_affine_tool - only one conversion at a time    

    Args:
        fixed_file_path (PosixPath): Path to fixed image
        moving_file_path (PosixPath): Path to moving image
        mat_file_path (PosixPath): Path to c3d matrix
        ras2fsl (bool, optional): convert C3D to FSL format. Defaults to True.
        ras2itk (bool, optional): convert C3D to ITK format. Defaults to False.
        fsl2ras (bool, optional): convert FSL to C3D format. Defaults to False.
        fsl2itk (bool, optional): convert FSL to ITK format. Defaults to False.
    """

    if type(fixed_file_path) is str:
        fixed_file_path = Path(fixed_file_path)
    elif type(moving_file_path) is str:
        moving_file_path = Path(moving_file_path)
    elif type(mat_file_path) is str:
        mat_file_path = Path(mat_file_path)

    out_prefix = misc.get_prefix(moving_file_path)

    pre_cmd = "c3d_affine_tool " + \
        "-ref " + fixed_file_path.as_posix() + \
        " " + \
        "-src " + moving_file_path.as_posix() + \
        " " + \
        mat_file_path.as_posix() + " "

    if ras2fsl is True:
        main_cmd = pre_cmd + \
            "-ras2fsl " + " -o " + \
            out_prefix.as_posix() + "_gdy2fsl.txt"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_prefix.as_posix() + "_gdy2fsl.txt"
    elif ras2itk is True:
        main_cmd = pre_cmd + \
            "-oitk " + " " + \
            out_prefix.as_posix() + "_gdy2itk.txt"
        outfile = out_prefix.as_posix() + "_gdy2itk.txt"
        misc.exec_shell(cmd=main_cmd)
    elif fsl2ras is True:
        main_cmd = pre_cmd + \
            "-fsl2ras " + " -o " + \
            out_prefix.as_posix() + "_fsl2gdy.txt"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_prefix.as_posix() + "_fsl2gdy.txt"
    elif fsl2itk is True:
        main_cmd = pre_cmd + \
            "-fsl2ras " + "-oitk " + \
            out_prefix.as_posix() + "_fsl2itk.txt"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_prefix.as_posix() + "_fsl2itk.txt"
    return outfile

def c3d_reslice(fixed_file_path, moving_file_path, mat_file_path, interp="Cubic", verbose=False):
    """Apply transformations using c3d

    Args:
        fixed_file_path (PosixPath): Path to the fixed image
        moving_file_path (PosixPath): Path to the moving image
        mat_file_path (PosixPath): Path to a c3d affine matrix
        interp (str, optional): "Cubic" or "Sinc". Defaults to "Cubic".
        verbose (bool, optional): [description]. Defaults to False.
    """

    out_file_prefix = misc.get_prefix(moving_file_path)
    out_file_path = Path(out_file_prefix.as_posix() + "_c3dWarped.nii.gz")

    main_cmd = "c3d " + \
        fixed_file_path.as_posix() + " " + \
        moving_file_path.as_posix() + " " + \
        "-interpolation " + interp + " " + \
        "-reslice-matrix " + mat_file_path.as_posix() + " " + \
        "-o " + out_file_path.as_posix()
    misc.exec_shell(cmd=main_cmd)

    return out_file_path


def greedy_register(fixed_file_path, moving_file_path, transform="rigid", cost="NMI", stages=2, init="identity", mask=None, verbose=False):
    """Register using greedy

    Args:
        fixed_file_path (PosixPath): Path to fixed image
        moving_file_path (PosixPath): Path to moving image
        transform (str, optional): "rigid" or "affine". Defaults to "rigid".
        cost (str, optional): "NMI" or "NCC". Defaults to "NMI".
        stages (int, optional): 4 --> 200x100x100x50 
                                3 --> 100x100x50 
                                2 --> 100x50
                                1 --> 50
                                Defaults to 2.
        init (PosixPath or str, optional): Path to initial matrix or "centres". Defaults to "identity".
        mask (PosixPath): Path to mask image. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to False.
    """
    num_cores = mp.cpu_count()
    out_file_prefix = misc.get_prefix(moving_file_path)

    if transform == "rigid":
        dof = str(6)
    elif transform == "affine":
        dof = str(12)

    if cost == "NCC":
        cost = str("NCC 4x4x4")

    if stages == 1:
        iters = str("50")
    elif stages == 2:
        iters = str("100x50")
    elif stages == 3:
        iters = str("100x100x50")
    elif stages == 4:
        iters = str("200x100x100x50")

    if init == "identity":
        initialisation = "-ia-identity"
    elif init == "centres":
        initialisation = "-ia-image-centers"
    else:
        initialisation = "-ia " + init.as_posix()

    if mask is None:
        main_cmd = "greedy -d 3 -a " + \
            initialisation + " " + \
            "-dof " + dof + " " + \
            "-threads " + str(num_cores) + " " + \
            "-m " + cost + " " + \
            "-n " + iters + " " + \
            "-i " + fixed_file_path.as_posix() + " " + moving_file_path.as_posix() + " " + \
            "-o " + out_file_prefix.as_posix() + "_gdy.mat"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_file_prefix.as_posix() + "_gdy.mat"
    else:
        main_cmd = "greedy -d 3 -a " + \
            initialisation + " " + \
            "-dof " + dof + " " + \
            "-threads " + str(num_cores) + " " + \
            "-m " + cost + " " + \
            "-n " + iters + " " + \
            "-i " + fixed_file_path.as_posix() + " " + moving_file_path.as_posix() + " " + \
            "-gm " + mask.as_posix() + " " + \
            "-o " + out_file_prefix.as_posix() + "_gdy.mat"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_file_prefix.as_posix() + "_gdy.mat"
    return outfile

def greedy_reslice(fixed_file_path, moving_file_path, mat_file_path, interp="label", invert=False):
    """Reslice using greedy 

    Args:
        fixed_file_path (PosixPath): [description]
        moving_file_path (PosixPath): [description]
        mat_file_path (PosixPath): [description]
        interp (str, optional): "nearest","linear" or "label". Defaults to "label".
        invert (bool, optional): invert affine matrix if True. Defaults to False.
    """
    num_cores = mp.cpu_count()
    out_file_prefix = misc.get_prefix(moving_file_path)
    
    if interp=="label":
        interp=str("LABEL 0.2vox")
    elif interp=="linear":
        interp="LINEAR"
    elif interp=="nearest":
        interp="NN"

    if invert is True:
        main_cmd = "greedy -d 3 -r " + mat_file_path.as_posix() + " ,-1 " + \
            "-ri " + interp + " " + \
            "-rf " + fixed_file_path.as_posix() + " " + \
            "-rm " + moving_file_path.as_posix() + " " + out_file_prefix.as_posix() + "_gdyWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)
        outfile=out_file_prefix.as_posix() + "_gdyWarped.nii.gz"
    else:
        main_cmd = "greedy -d 3 -r " + mat_file_path.as_posix() + " " + \
            "-ri " + interp + " " + \
            "-rf " + fixed_file_path.as_posix() + " " + \
            "-rm " + moving_file_path.as_posix() + " " + out_file_prefix.as_posix() + "_gdyWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)
        outfile=out_file_prefix.as_posix() + "_gdyWarped.nii.gz"
    return outfile