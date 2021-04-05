#!/bin/env python
"""
Python scripts for FSL command-line utilities
"""
import os
import sys
from . import misc
import shutil
import subprocess
import numpy as np
from pathlib import Path
import multiprocessing as mp


def fsl_split(input_file_path, out_prefix=None, verbose=False):
    """Wraps fslsplit
    Args:
        input_file_path (PosixPath): Path to the input file
    Returns:
        PosixPath: path to the output files without extension
    """

    misc.make_dir(input_file_path)

    if out_prefix is None:
        out_prefix = misc.get_prefix(input_file_path)
        out_prefix = out_prefix / "vol_"

    main_cmd = "fslsplit " + \
        input_file_path.as_posix() + " " + \
        out_prefix.as_posix() + " " + \
        "-t"

    misc.exec_shell(cmd=main_cmd)

    return out_prefix


def fsl_merge(input_file_list, out_file_path, tr=1.0, verbose=False):
    """Wraps fslmerge
    Args:
        input_file_list (list): List of input files
        out_file_path (PosixPath): Path to the output file
        tr (float, optional): TR of the data. Defaults to 1.0.
    """
    if input_file_list is list:
        input_file_string = misc.list_to_string(in_list=input_file_list)
    else:
        input_file_string = input_file_list

    if out_file_path is str:
        out_file_path = Path(out_file_path)

    main_cmd = "fslmerge -tr " + \
        out_file_path.as_posix() + " " + \
        input_file_string + " " + \
        str(tr)

    misc.exec_shell(cmd=main_cmd)


def fsl_register(fixed_file_path, moving_file_path, transform="rigid", cost="NMI", init=None, run=True):
    """Rigid register using flirt 

    Args:
        fixed_file_path (PosixPath): Path to fixed file.
        moving_file_path (PosixPath): Path to moving file.
        transform (str, optional): "rigid" or "affine". Defaults to "rigid".
        cost (str, optional):   "NMI" --> "normmi"
                                "NC"  --> "normcorr"
                                "CR"  --> "corrratio"
                                Defaults to "NMI".
        init (PosixPath, optional): Path to initial matrix. Defaults to None.
        run (str, optional): Set to False to only export command. Defaults to True.
    """
    out_file_prefix = misc.get_prefix(moving_file_path)

    if transform == "rigid":
        dof = 6
    elif transform == "affine":
        dof = 12

    if cost == "NMI":
        cost = str("normmi")
    elif cost == "NC":
        cost = str("normcorr")
    elif cost == "CR":
        cost = str("corratio")

    if init is None:
        initialisation = " "
    else:
        initialisation = "-init " + init.as_posix() + " "

    main_cmd = "flirt " + \
        "-interp spline" + " " + \
        "-bins 512" + " " + \
        "-noresample " + \
        "-noclamp " + \
        "-noresampblur " + \
        "-nosearch " + \
        "-cost " + cost + " " + \
        "-dof " + str(dof) + " " + \
        initialisation + \
        "-ref" + " " + \
        fixed_file_path.as_posix() + " " + \
        "-in" + " " + \
        moving_file_path.as_posix() + " " + \
        "-out" + " " + \
        out_file_prefix.as_posix() + "_fslWarped.nii.gz" + " " + \
        "-omat" + " " + \
        out_file_prefix.as_posix() + "_fsl.mat"
    if run is True:
        misc.exec_shell(cmd=main_cmd)
        out_file = out_file_prefix.as_posix() + "_fslWarped.nii.gz"
        return out_file
    else:
        return main_cmd


def fsl_realign(input_file_list, fixed_file_path=None):
    """Run flirt in parallel for 3D rigid realignment

    Args:
        input_file_list (PosixPath): List of files for realignment
        fixed_file_path (PosixPath, optional): Path to a reference file. Defaults to None.
    """
    realign_cmd_list = list()

    if fixed_file_path is None:
        for vol_num in range(1, len(input_file_list)):
            cmd = fsl_register(fixed_file_path=input_file_list[0],
                               moving_file_path=input_file_list[vol_num], run=False)
            realign_cmd_list.append(cmd)
        misc.exec_parashell(cmd_list=realign_cmd_list)
    else:
        for vol_num in range(0, len(input_file_list)):
            cmd = fsl_register(fixed_file_path=fixed_file_path,
                               moving_file_path=input_file_list[vol_num], run=False)
            realign_cmd_list.append(cmd)
        misc.exec_parashell(cmd_list=realign_cmd_list)


def fsl_reslice(fixed_file_path, moving_file_path, mat_file_path, interp="spline"):
    """Apply estimated transformations using flirt

    Args:
        fixed_file_path (PosixPath): Path to fixed file.
        moving_file_path (PosixPath): Path to moving file.
        mat_file_path (PosixPath): Path to affine matrix.
        interp (str, optional): "spline" or "sinc". Defaults to "spline".
    """
    out_file_prefix = misc.get_prefix(moving_file_path)

    main_cmd = "flirt " + \
        "-applyxfm" + " " + \
        "-interp " + interp + " " + \
        "-init " + mat_file_path.as_posix() + " " + \
        "-ref " + fixed_file_path.as_posix() + " " + \
        "-in " + moving_file_path.as_posix() + " " + \
        "-out " + out_file_prefix.as_posix() + "_fslWarped.nii.gz"

    misc.exec_shell(cmd=main_cmd)
    outfile = out_file_prefix.as_posix() + "_fslWarped.nii.gz"
    return outfile


def fsl_invertmat(mat_file_path):
    """Invert flirt matrix

    Args:
        mat_file_path (PosixPath): Path to affine matrix.
    """
    out_file_prefix = misc.get_prefix(mat_file_path)

    main_cmd = "convert_xfm " + \
        "-omat " + \
        out_file_prefix.as_posix() + "_inverse.mat" + " " + \
        "-inverse " + \
        mat_file_path.as_posix()
    misc.exec_shell(cmd=main_cmd)
    return out_file_prefix.as_posix() + "_inverse.mat"


def fsl_concatmat(mat_A2B_path, mat_B2C_path, mat_A2C_out):
    """Concat flirt matrix to get mat_A2C_out

    Args:
        mat_A2B_path (PosixPath): Path to affine matrix.
        mat_B2C_path (PosixPath): Path to affine matrix.
        mat_A2C_out (PosixPath): Path to affine matrix.        
    """
    main_cmd = "convert_xfm " + \
        "-omat " + \
        mat_A2C_out.as_posix() + " " + \
        "-concat " + \
        mat_B2C_path.as_posix() + " " + \
        mat_A2B_path.as_posix()
    misc.exec_shell(cmd=main_cmd)


def fsl_avscale(mat_file_path):
    """Wrap avscale. Obtain rotation and translation info from FSL mat file.

    Args:
        mat_file_path (PosixPath): [description]
    """
    main_cmd = "avscale --allparams " + mat_file_path.as_posix()

    mat_file_output = subprocess.check_output(main_cmd, shell=True)

    return mat_file_output


def fsl_Tmean(input_file_path):
    """Calculate temporal mean
    Args:
        input_file_path (PosixPath): Path to the input file
    """
    out = misc.get_prefix(input_file_path)
    out_string = out.as_posix()+"_fslTmean.nii.gz"

    main_cmd = "fslmaths " + \
        input_file_path.as_posix() + " " + \
        "-Tmean" + " " + \
        out_string

    misc.exec_shell(cmd=main_cmd)

    return Path(out_string)


def fsl_Tstd(input_file_path):
    """Calculate temporal standard deviation
    Args:
        input_file_path (PosixPath): Path to the input file
    """
    out = misc.get_prefix(input_file_path)
    out_string = out.as_posix()+"_fslTstd.nii.gz"

    main_cmd = "fslmaths " + \
        input_file_path.as_posix() + " " + \
        "-Tstd" + " " + \
        out_string

    misc.exec_shell(cmd=main_cmd)

    return Path(out_string)


def fsl_Tsnr(input_file_path):
    """Calculate temporal SNR after calculating Tmean and Tstd
    Args:
        input_file_path (PosixPath): Path to the input file
    """
    out = misc.get_prefix(input_file_path)
    out_string = out.as_posix()+"_fslTsnr.nii.gz"

    fsl_Tmean(input_file_path=input_file_path)
    fsl_Tstd(input_file_path=input_file_path)

    main_cmd = "fslmaths" + " " + \
        out.as_posix() + "_fslTmean.nii.gz" + " " + \
        "-div" + " " + \
        out.as_posix() + "_fslTstd.nii.gz" + " " + \
        out_string

    misc.exec_shell(cmd=main_cmd)

    return Path(out_string)


def fsl_bet(input_file_path, frac=0.25):
    """Wrap bet, basic settings

    Args:
        input_file_path (PosixPath): Path to the input file
        frac (float, optional): Fractional threshold (larger mask < 0.5 < smaller mask). Defaults to 0.25.
    """
    if type(input_file_path) is str:
        src = Path(input_file_path)
    out = misc.get_prefix(input_file_path)

    main_cmd = "bet " + \
        input_file_path.as_posix() + " " + \
        out.as_posix() + "_fslBet" + " " + \
        "-f 0.25 -m"

    misc.exec_shell(cmd=main_cmd)


def fsl_topup(input_file_path, topup_parameters, topup_config=None, verbose=False):
    """Wrap topup

    Args:
        input_file_path (PosixPath): Path to 4D file to estimate topup
        topup_parameters (PosixPath): Path to topup parameters file
        topup_config (PosixPath, optional): Path to topup configuration. Defaults to topup default.
    """
    if type(input_file_path) is str:
        input_file_path = Path(input_file_path)
    out = misc.get_prefix(input_file_path)

    if topup_config is not None:
        main_cmd = "topup " + \
            "--config=" + topup_config.as_posix() + " " + \
            "--out=" + out.as_posix() + "_fslTopup" + " " + \
            "--imain=" + input_file_path.as_posix() + " " + \
            "--datain=" + topup_parameters.as_posix() + " " + \
            "--fout=" + out.as_posix() + "_fslTopup_fieldmap_in_Hz.nii.gz" + " " + \
            "--iout=" + out.as_posix() + "_fslTopup_unWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)
    else:
        main_cmd = "topup " + \
            "--out=" + out.as_posix() + "_fslTopup" + " " + \
            "--imain=" + input_file_path.as_posix() + " " + \
            "--datain=" + topup_parameters.as_posix() + " " + \
            "--fout=" + out.as_posix() + "_fslTopup_fieldmap_in_Hz.nii.gz" + " " + \
            "--iout=" + out.as_posix() + "_fslTopup_unWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)

    # Get fieldmap magnitude
    main_cmd = "fslmaths " + \
        out.as_posix() + "_fslTopup_unWarped.nii.gz" + " " + \
        "-Tmean " + \
        out.as_posix() + "_fslTopup_fieldmap_mag.nii.gz"

    misc.exec_shell(cmd=main_cmd)

    # Bet fieldmap magnitude
    fsl_bet(input_file_path=Path(
        out.as_posix() + "_fslTopup_fieldmap_mag.nii.gz"))

    # Rescale fieldmap from Hz to radians/sec
    main_cmd = "fslmaths " + \
        out.as_posix() + "_fslTopup_fieldmap_in_Hz.nii.gz" + " " + \
        "-mul 6.28319" + " " + \
        out.as_posix() + "_fslTopup_fieldmap_in_radPersec.nii.gz"

    misc.exec_shell(cmd=main_cmd)

    output_dict = {"fieldmap_Hz": out.as_posix() + "_fslTopup_fieldmap_in_Hz.nii.gz",
                   "fieldmap_radPersec": out.as_posix() + "_fslTopup_fieldmap_in_radPersec.nii.gz",
                   "fieldmap_mag": out.as_posix() + "_fslTopup_fieldmap_mag.nii.gz",
                   "fieldmap_mag_brain": out.as_posix() + "_fslTopup_fieldmap_mag_fslBet.nii.gz",
                   "fieldmap_mag_brain_mask": out.as_posix() + "_fslTopup_fieldmap_mag_fslBet_mask.nii.gz"}
    return output_dict


def fsl_applytopup(input_file_path, topup_parameters, topup_output_prefix, in_index=1, verbose=False):
    """Wrap applytopup

    Args:
        input_file_path (PosixPath): Path to 4D file to estimate topup
        topup_parameters (PosixPath): Path to topup parameters file
        topup_output_prefix (str): Path to topup output prefix
        in_index(int, optional): Input file index for unwarping
    """
    if type(input_file_path) is str:
        input_file_path = Path(input_file_path)
    out = misc.get_prefix(input_file_path)

    if type(topup_output_prefix) is not str:
        topup_output_prefix = topup_output_prefix.as_posix()

    main_cmd = "applytopup " + \
        "--method=jac" + " " + \
        "--out=" + out.as_posix() + "_fslDiCorr" + " " + \
        "--imain=" + input_file_path.as_posix() + " " + \
        "--datain=" + topup_parameters.as_posix() + " " + \
        "--topup=" + topup_output_prefix + " " + \
        "--inindex=" + str(in_index)

    misc.exec_shell(cmd=main_cmd)

    # Calc temporal mean
    fsl_Tmean(input_file_path=Path(out.as_posix() + "_fslDiCorr.nii.gz"))


def fsl_surrDiff(input_file_path, order="tc", tr=1.0, add_mean=True, bold=False, verbose=False):
    """Do surround subtraction the FSL way command

    Args:
        input_file_path (PosixPath): Path to input file.
        order (str, optional): "tc" or "ct". Defaults to "tc".
        tr (float, optional): Volume TR. Defaults to 1.0.
        add_mean (bool, optional): Rescale. Defaults to True.
        bold (bool, optional): Also calculate BOLD timeseries. Defaults to False.
        verbose (bool, optional): Print verbose output. Defaults to False.
    """
    if type(input_file_path) is str:
        input_file_path = Path(input_file_path)
    out_prefix = misc.get_prefix(input_file_path)

    if order == "ct":
        invert = "-mul -1.0"

    misc.make_dir(input_file_path=input_file_path)

    fsl_split(input_file_path=input_file_path)

    even_filelist = list(out_prefix.glob("vol_*[02468].nii.gz"))
    even_filelist.sort()
    even_filelist = misc.list_to_string(even_filelist)
    odd_filelist = list(out_prefix.glob("vol_*[13579].nii.gz"))
    odd_filelist.sort()
    odd_filelist = misc.list_to_string(odd_filelist)

    # For when the timeseries has odd number of TRs
    if len(even_filelist) != len(odd_filelist):
        even_filelist = even_filelist[:-1]
        fsl_merge(input_file_list=even_filelist,
                  out_file_path=Path(out_prefix.as_posix() + "/even.nii.gz"))
    else:
        fsl_merge(input_file_list=even_filelist,
                  out_file_path=Path(out_prefix.as_posix() + "/even.nii.gz"))

    fsl_merge(input_file_list=odd_filelist,
              out_file_path=Path(out_prefix.as_posix() + "/odd.nii.gz"))

    def temp_shift(infile, outfile, shift):
        """Temporal interpolation using slicetimer

        Args:
            infile (str): Path to input file.
            outfile (str): Path to output file.
            shift (float): shift factor.
        """
        main_cmd = "slicetimer " + \
            "-i " + infile + " " + \
            "-o " + outfile + " " + \
            "--tglobal=" + str(shift)

        misc.exec_shell(cmd=main_cmd)

    temp_shift(infile=out_prefix.as_posix() + "/odd.nii.gz",
               outfile=out_prefix.as_posix() + "/odd_earlier.nii.gz",
               shift=0.25)
    temp_shift(infile=out_prefix.as_posix() + "/odd.nii.gz",
               outfile=out_prefix.as_posix() + "/odd_later.nii.gz",
               shift=-0.25)
    temp_shift(infile=out_prefix.as_posix() + "/even.nii.gz",
               outfile=out_prefix.as_posix() + "/even_earlier.nii.gz",
               shift=0.25)
    temp_shift(infile=out_prefix.as_posix() + "/even.nii.gz",
               outfile=out_prefix.as_posix() + "/even_later.nii.gz",
               shift=0.75)

    main_cmd = "fslmaths " + \
        out_prefix.as_posix() + "/odd_later.nii.gz" + " " + \
        "-sub " + \
        out_prefix.as_posix() + "/even_earlier.nii.gz" + " " + \
        out_prefix.as_posix() + "/even_to_odd.nii.gz"

    misc.exec_shell(cmd=main_cmd)

    main_cmd = "fslmaths " + \
        out_prefix.as_posix() + "/odd_earlier.nii.gz" + " " + \
        "-sub " + \
        out_prefix.as_posix() + "/even_later.nii.gz" + " " + \
        out_prefix.as_posix() + "/odd_to_even.nii.gz"

    misc.exec_shell(cmd=main_cmd)

    surr_diff_out_path = Path(out_prefix.as_posix() + "/surrDiff")

    misc.make_dir(input_file_path=surr_diff_out_path)

    even_sub_path = Path(surr_diff_out_path.as_posix() + "/e_")
    fsl_split(input_file_path=Path(out_prefix.as_posix() + "/even_to_odd.nii.gz"),
              out_prefix=even_sub_path)
    even_sub_filelist = list(surr_diff_out_path.glob("e_*.nii.gz"))
    even_sub_filelist.sort()

    odd_sub_path = Path(surr_diff_out_path.as_posix() + "/o_")
    fsl_split(input_file_path=Path(out_prefix.as_posix() + "/odd_to_even.nii.gz"),
              out_prefix=odd_sub_path)

    odd_sub_filelist = list(surr_diff_out_path.glob("o_*.nii.gz"))
    odd_sub_filelist.sort()

    subtracted_fulllist = even_sub_filelist + odd_sub_filelist
    subtracted_fulllist[::2] = even_sub_filelist
    subtracted_fulllist[1::2] = odd_sub_filelist

    subtracted_fulllist_merge = misc.list_to_string(
        in_list=subtracted_fulllist)
    fsl_merge(input_file_list=subtracted_fulllist_merge,
              out_file_path=Path(out_prefix.as_posix() + "_surrDiff.nii.gz"),
              tr=tr)
    fsl_Tmean(input_file_path=Path(out_prefix.as_posix() + "_surrDiff.nii.gz"))

    if add_mean is True:
        mean_out = fsl_Tmean(input_file_path)

        main_cmd = "fslmaths " + \
            out_prefix.as_posix() + "_surrDiff.nii.gz" + " " + \
            "-add " + \
            out_prefix.as_posix() + "_fslTmean.nii.gz" + " " + \
            out_prefix.as_posix() + "_surrDiff_plusMean.nii.gz"
        misc.exec_shell(cmd=main_cmd)

    if bold is True:
        main_cmd = "fslmaths " + \
            out_prefix.as_posix() + "/odd_later.nii.gz" + " " + \
            "-add " + \
            out_prefix.as_posix() + "/even_earlier.nii.gz" + " " + \
            "-div 2 " + \
            out_prefix.as_posix() + "/even_to_odd.nii.gz"
        misc.exec_shell(cmd=main_cmd)

        main_cmd = "fslmaths " + \
            out_prefix.as_posix() + "/odd_earlier.nii.gz" + " " + \
            "-add " + \
            out_prefix.as_posix() + "/even_later.nii.gz" + " " + \
            "-div 2 " + \
            out_prefix.as_posix() + "/odd_to_even.nii.gz"
        misc.exec_shell(cmd=main_cmd)

        surr_avg_out_path = Path(out_prefix.as_posix() + "/surrAvg")
        misc.make_dir(input_file_path=surr_avg_out_path)

        even_sub_path = Path(surr_avg_out_path.as_posix() + "/e_")
        fsl_split(input_file_path=Path(out_prefix.as_posix() + "/even_to_odd.nii.gz"),
                  out_prefix=even_sub_path)
        even_sub_filelist = list(surr_avg_out_path.glob("e_*.nii.gz"))
        even_sub_filelist.sort()

        odd_sub_path = Path(surr_avg_out_path.as_posix() + "/o_")
        fsl_split(input_file_path=Path(out_prefix.as_posix() + "/odd_to_even.nii.gz"),
                  out_prefix=odd_sub_path)
        odd_sub_filelist = list(surr_avg_out_path.glob("o_*.nii.gz"))
        odd_sub_filelist.sort()

        averaged_fulllist = even_sub_filelist + odd_sub_filelist
        averaged_fulllist[::2] = even_sub_filelist
        averaged_fulllist[1::2] = odd_sub_filelist

        averaged_fulllist_merge = misc.list_to_string(
            in_list=averaged_fulllist)
        fsl_merge(input_file_list=averaged_fulllist_merge,
                  out_file_path=Path(
                      out_prefix.as_posix() + "_surrAvg.nii.gz"),
                  tr=tr)
        fsl_Tmean(input_file_path=Path(
            out_prefix.as_posix() + "_surrAvg.nii.gz"))
        misc.remove_dir(surr_avg_out_path)

    misc.remove_files(src=out_prefix,
                      pattern="vol*.nii.gz")
    misc.remove_files(src=out_prefix,
                      pattern="odd*.nii.gz")
    misc.remove_files(src=out_prefix,
                      pattern="even*.nii.gz")
    misc.remove_dir(surr_diff_out_path)
    misc.remove_dir(out_prefix)


def fsl_bbr(fixed_file_path, moving_file_path, wm_seg_path, interp="spline", init_file_path=None):
    """Register using BBR

    Args:
        fixed_file_path ([type]): Path to fixed file.
        moving_file_path ([type]): Path to moving file.
        wm_seg_path ([type]): Path to WM Seg file.
        interp (str, optional): "spline" or "sinc". Defaults to "spline".
        init_file_path (PosixPath, optional): Path to init mat file. Defaults to None.
    """

    out_file_prefix = misc.get_prefix(moving_file_path)
    if init_file_path is not None:
        main_cmd = "flirt -nosearch -noresample -noresampblur " + \
            "-cost bbr " + \
            "-interp " + interp + " " + \
            "-init " + init_file_path.as_posix() + " " + \
            "-ref " + fixed_file_path.as_posix() + " " + \
            "-wmseg " + wm_seg_path.as_posix() + " " + \
            "-in " + moving_file_path.as_posix() + " " + \
            "-omat " + out_file_prefix.as_posix() + "_fslBBR.mat" + " " + \
            "-out " + out_file_prefix.as_posix() + "_fslbbWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_file_prefix.as_posix() + "_fslBBR.mat"
    else:
        main_cmd = "flirt -nosearch -noresample -noresampblur " + \
            "-schedule ${FSLDIR}/etc/flirtsch/bbr.sch -cost bbr " + \
            "-interp " + interp + " " + \
            "-ref " + fixed_file_path.as_posix() + " " + \
            "-wmseg " + wm_seg_path.as_posix() + " " + \
            "-in " + moving_file_path.as_posix() + " " + \
            "-omat " + out_file_prefix.as_posix() + "_fslBBR.mat" + " " + \
            "-out " + out_file_prefix.as_posix() + "_fslbbWarped.nii.gz"
        misc.exec_shell(cmd=main_cmd)
        outfile = out_file_prefix.as_posix() + "_fslBBR.mat"
    return outfile


def fsl_despike(input_file_path, mask_file_path):
    """Despike fieldmap

    Args:
        input_file_path (PosixPath): Path to input file.
        mask_file_path (PosixPath): Path to mask file.
    """
    if type(input_file_path) is str:
        input_file_path = Path(input_file_path)
    if type(mask_file_path) is str:
        mask_file_path = Path(mask_file_path)

    out_prefix = misc.get_prefix(input_file_path)
    out_mask_prefix = misc.get_prefix(mask_file_path)

    main_cmd = "fslmaths " + \
        mask_file_path.as_posix() + " " + \
        "-ero " + \
        out_mask_prefix.as_posix() + "_ero.nii.gz"

    misc.exec_shell(cmd=main_cmd)

    main_cmd = "fugue " + \
        "--loadfmap=" + input_file_path.as_posix() + " " + \
        "--savefmap=" + out_prefix.as_posix() + "_filt.nii.gz" + " " + \
        "--mask=" + out_mask_prefix.as_posix() + "_ero.nii.gz" + " " + \
        "--despike" + " " + "--despikethreshold=2.1"
    misc.exec_shell(cmd=main_cmd)

    main_cmd = "fslmaths " + \
        input_file_path.as_posix() + " " + \
        "-sub " + \
        out_prefix.as_posix() + "_filt.nii.gz" + " " + \
        "-mas " + mask_file_path.as_posix() + " " + \
        "-add " + out_prefix.as_posix() + "_filt.nii.gz" + " " + \
        out_prefix.as_posix() + "_despiked.nii.gz"
    misc.exec_shell(cmd=main_cmd)

    file_rm = Path(out_prefix.as_posix() + "_filt.nii.gz")
    file_rm.unlink()
    file_rm = Path(out_mask_prefix.as_posix() + "_ero.nii.gz")
    file_rm.unlink()


def fsl_eddy(input_file_path, bvecs_file, bvals_file, slspec_file, mask_file, topup_parameters, index_file, topup_out, mt="gpu"):
    """Wrap eddy_cuda9.1

    Args:
        input_file_path (PosixPath): Path to 4D file
        bvecs_file (PosixPath): Path to bvecs file
        bvals_file (PosixPath): Path to bvals file
        slspec_file (PosixPath): Path to slspec file
        mask_file (PosixPath): Path to brain mask file
        topup_parameters (PosixPath): Path to topup parameters file
        index_file (PosixPath): Path to index file.
        topup_out (PosixPath, optional): Path to topup output prefix
        mt (str, optional): controls multi-threading. Defaults to "gpu".
                            can be: cpu, None

    """
    if type(input_file_path) is str:
        input_file_path = Path(input_file_path)
    out_prefix = misc.get_prefix(input_file_path)

    if mt == "gpu":
        eddy_type = "eddy_cuda9.1"
        main_cmd = eddy_type + " " + \
            "--imain=" + input_file_path.as_posix() + " " + \
            "--bvecs=" + bvecs_file.as_posix() + " " + \
            "--bvals=" + bvals_file.as_posix() + " " + \
            "--slspec=" + slspec_file.as_posix() + " " + \
            "--mask=" + mask_file.as_posix() + " " + \
            "--acqp=" + topup_parameters.as_posix() + " " + \
            "--index=" + index_file.as_posix() + " " + \
            "--topup=" + topup_out.as_posix() + " " + \
            "--out=" + out_prefix.as_posix() + "_fslEddy"
    elif mt == "cpu":
        # CPU versions of eddy do not have slspec functionality
        eddy_type = "eddy_openmp"
        main_cmd = eddy_type + " " + \
            "--imain=" + input_file_path.as_posix() + " " + \
            "--bvecs=" + bvecs_file.as_posix() + " " + \
            "--bvals=" + bvals_file.as_posix() + " " + \
            "--mask=" + mask_file.as_posix() + " " + \
            "--acqp=" + topup_parameters.as_posix() + " " + \
            "--index=" + index_file.as_posix() + " " + \
            "--topup=" + topup_out.as_posix() + " " + \
            "--out=" + out_prefix.as_posix() + "_fslEddy"
    else:
        # CPU versions of eddy do not have slspec functionality
        eddy_type = "eddy"
        main_cmd = eddy_type + " " + \
            "--imain=" + input_file_path.as_posix() + " " + \
            "--bvecs=" + bvecs_file.as_posix() + " " + \
            "--bvals=" + bvals_file.as_posix() + " " + \
            "--mask=" + mask_file.as_posix() + " " + \
            "--acqp=" + topup_parameters.as_posix() + " " + \
            "--index=" + index_file.as_posix() + " " + \
            "--topup=" + topup_out.as_posix() + " " + \
            "--out=" + out_prefix.as_posix() + "_fslEddy"

    misc.exec_shell(cmd=main_cmd)
