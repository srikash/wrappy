#!/bin/env python
"""
Python scripts for MRTRIX 3 command-line utilities
"""
import os
import sys
from . import misc
from pathlib import Path

def gibbs_unring(input_file_path):
    """ Remove Gibbs ringing artefacts from data

    Args:
        input_file_path (PosixPath): Path to input file
    """
    out_file_prefix = misc.get_prefix(input_file_path)

    main_cmd = "mrdegibbs" + " " + \
        input_file_path + " " + \
        out_file_prefix.as_posix() + "deG.nii.gz"

def denoise_diffusion(input_file_path):
    """MP-PCA denoising for DWI

    Args:
        input_file_path (PosixPath): Path to input file
    """

    out_file_prefix = misc.get_prefix(input_file_path)

    main_cmd = "dwidenoise" + " " + \
        input_file_path.as_posix() + " " + \
        out_file_prefix.as_posix() + "deN.nii.gz"

    misc.exec_shell(cmd=main_cmd)
