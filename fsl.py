#!/bin/env python
"""
Python scripts for FSL command-line utilities
"""
import os
import sys
import misc
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
        out_prefix + " " + \
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
    input_file_string = misc.list_to_string(in_list=input_file_list)
    
    main_cmd = "fslmerge -tr " + \
        out_file_path.as_posix() + " " + \
        input_file_string + " " + \
        str(tr)

    misc.exec_shell(cmd=main_cmd)


def fsl_register(fixed_file_path, moving_file_path, transform="rigid", cost="NMI",init=None):
    """[summary]

    Args:
        fixed_file_path ([type]): [description]
        moving_file_path ([type]): [description]
        transform (str, optional): [description]. Defaults to "rigid".
        cost (str, optional): [description]. Defaults to "NMI".
        init ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
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
        
    misc.exec_shell(cmd=main_cmd)


def fsl_realign(input_file_list, ref_file_path=None):
    """[summary]

    Args:
        input_file_list ([type]): [description]
        ref_file_path ([type], optional): [description]. Defaults to None.
    """
    
    realign_cmd_list = list()
    if ref_file_path is None:
        for vol_num in range(1, len(input_file_list)):
            cmd = fsl_register(fixed_file_path=input_file_list[0], 
                               moving_file_path=input_file_list[vol_num])
            realign_cmd_list.append(cmd)
        misc.exec_parashell(moco_cmd_list)
    else:
        for vol_num in range(len(input_file_list)):
            cmd = fsl_register(fixed_file_path=ref_file_path, 
                        moving_file_path=input_file_list[vol_num])
            realign_cmd_list.append(cmd)
        misc.exec_parashell(realign_cmd_list)

def fsl_avscale(mat_file_path):
    """[summary]

    Args:
        mat_file_path ([type]): [description]
    """
    main_cmd = "avscale --allparams " + \
        str(mat_file_path)
    mat_file_output = subprocess.check_output(main_cmd, shell=True)
    return mat_file_output

def fsl_Tmean(src, verbose=False):
    """Wrapper for fslmaths command
    Args:
        src (PosixPath): path to the input files
    Returns:
        PosixPath: path to the output files without extension
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Calculating Tmean")
        print("%---------------------------------------%")
    if type(src) is str:
        src = Path(src)
    out = get_prefix(src)
    out_string = str(out)+"_Tmean.nii.gz"

    main_cmd = "fslmaths " + \
        str(src) + \
        " " + \
        "-Tmean " + \
        out_string

    run_cmd = main_cmd
    if verbose is True:
        print("++++ Running fslmaths ++++")
        print(run_cmd)
        print("%---------------------------------------%")
    subprocess.run(run_cmd, shell=True)
    return Path(out_string)


def calc_Tstd(src, verbose=False):
    """Wrapper for fslmaths command
    Args:
        src (PosixPath): path to the input files
    Returns:
        PosixPath: path to the output files without extension
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Calculating Tstd")
        print("%---------------------------------------%")
    if type(src) is str:
        src = Path(src)
    out = get_prefix(src)
    out_string = str(out)+"_Tstd.nii.gz"

    main_cmd = "fslmaths " + \
        str(src) + \
        " " + \
        "-Tstd " + \
        out_string

    run_cmd = main_cmd
    if verbose is True:
        print("++++ Running fslmaths ++++")
        print(run_cmd)
        print("%---------------------------------------%")
    subprocess.run(run_cmd, shell=True)
    return Path(out_string)


def calc_Tsnr(src, verbose=False):
    """Wrapper for fslmaths command
    Args:
        src (PosixPath): path to the input files
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Calculating Tmean")
        print("%---------------------------------------%")
    if type(src) is str:
        src = Path(src)

    out = get_prefix(src)
    calc_Tmean(src)
    calc_Tstd(src)
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Calculating calc_Tsnr")
        print("%---------------------------------------%")
    main_cmd = "fslmaths " + \
        str(out)+"_Tmean.nii.gz" + \
        " " + \
        "-div " + \
        str(out)+"_Tstd.nii.gz" + \
        " " + \
        str(out)+"_Tsnr.nii.gz"

    run_cmd = main_cmd
    if verbose is True:
        print("++++ Running fslmaths ++++")
        print(run_cmd)
        print("%---------------------------------------%")
    subprocess.run(run_cmd, shell=True)

def brainmask(src):
    """[summary]

    Args:
        src ([type]): [description]
    """
    if type(src) is str:
        src = Path(src)
    out = get_prefix(src)
    
    main_cmd = "bet " + \
        str(src) + \
        " " + \
        str(out) + "_brain" + \
        " "
    post_cmd = "-f 0.25 " + \
        "-m "

    run_cmd = main_cmd+post_cmd
    subprocess.run(run_cmd, shell=True)

def topup(imain, datain, config, verbose=False):
    """[summary]

    Args:
        imain ([type]): [description]
        datain ([type]): [description]
        config ([type]): [description]
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Distortion-correction using TOPUP")
        print("%---------------------------------------%")
    if type(imain) is str:
        imain = Path(imain)
    out = get_prefix(imain)

    if config is not None:
        pre_cmd = "topup " + \
            "--config=" + \
            str(config) + \
            " " + \
            "--out=" + \
            str(out)+"_topup "
    else:
        pre_cmd = "topup " + \
            "--out=" + \
            str(out)+"_topup "

    main_cmd = "--imain=" + \
        str(imain) + \
        " " + \
        "--datain=" + \
        str(datain) + \
        " " + \
        "--fout=" + \
        str(out)+"_topup_fmap_in_Hz.nii.gz " + \
        "--iout=" + \
        str(out)+"_topup_unwarped.nii.gz"

    run_cmd = pre_cmd+main_cmd
    if verbose is True:
        print("++++ Running TOPUP ++++")
        print(run_cmd)
    subprocess.run(run_cmd, shell=True)

    cmd_fmap_mag = "fslmaths " + \
        str(out)+"_topup_unwarped.nii.gz " + \
        "-Tmean " + \
        str(out)+"_topup_fmap_mag.nii.gz"
    if verbose is True:
        print("++++ Calculating fieldmap mag ++++")
        print(cmd_fmap_mag)
    subprocess.run(cmd_fmap_mag, shell=True)
    
    brainmask(src=Path(str(out)+"_topup_fmap_mag.nii.gz"))
    
    cmd_fmap_rad_per_s = "fslmaths " + \
    str(out)+"_topup_fmap_in_Hz.nii.gz " + \
    "-mul 6.28319" + \
    str(out)+"_topup_fmap_in_rad_per_s.nii.gz"
    
    if verbose is True:
        print("++++ Converting fieldmap to rad/s ++++")
        print(cmd_fmap_rad_per_s)
    subprocess.run(cmd_fmap_rad_per_s, shell=True)    

    cmd_fmap_rescale = "fslmaths " + \
        str(out)+"_topup_fmap_in_Hz.nii.gz " + \
        "-mul 6.28319 " + \
        str(out)+"_topup_fmap_in_rad_per_s.nii.gz"

    if verbose is True:
        print("++++ Rescaling fieldmap to rad/s ++++")
        print(cmd_fmap_rescale)
        print("%---------------------------------------%")
    subprocess.run(cmd_fmap_rescale, shell=True)

    output_dict={"fmap_hz": str(out)+"_topup_fmap_in_Hz.nii.gz",
                "fmap_rad_per_s": str(out)+"_topup_fmap_in_rad_per_s.nii.gz",
                "fmap_mag": str(out)+"_topup_fmap_mag.nii.gz",
                "fmap_mag_brain": str(out)+"_topup_fmap_mag_brain.nii.gz",
                "fmap_mag_brain_mask": str(out)+"_topup_fmap_mag_brain_mask.nii.gz"}
    return output_dict

def calc_surrDiff(src, first="t", tr=1.0, remean=True, surr_avg=False, verbose=False):
    """Do surround subtraction the FSL way command
    Args:
        src (PosixPath): path to the input files
    """
    if verbose is True:
        print("%---------------------------------------%")
        print("++ Calculating Surround Difference")
        print("%---------------------------------------%")
    if type(src) is str:
        src = Path(src)
    src_prefix = get_prefix(src)

    if first is "c":
        print("++++ First image is control ++++")
        invert = "-mul -1.0"

    make_dir(src)
    split4d(src)
    src_prefix = get_prefix(src)

    even_filelist = list(src_prefix.glob("vol_*[02468].nii.gz"))
    even_filelist.sort()
    even_filelist = list_to_string(even_filelist)
    odd_filelist = list(src_prefix.glob("vol_*[13579].nii.gz"))
    odd_filelist.sort()
    odd_filelist = list_to_string(odd_filelist)

    if len(even_filelist) != len(odd_filelist):
        print("++++ Ignoring last timepoint ++++")
        even_filelist = even_filelist[:-1]
        print(" ++ Merging to 4D file")
        merge3d(src=even_filelist,
                out=str(src_prefix)+"/even.nii.gz",
                tr=1.0)
    else:
        print(" ++ Merging to 4D file")
        merge3d(src=even_filelist,
                out=str(src_prefix)+"/even.nii.gz",
                tr=1.0)
    merge3d(src=odd_filelist,
            out=str(src_prefix)+"/odd.nii.gz",
            tr=1.0)

    def temp_shift(infile, outfile, shift):
        main_cmd = "slicetimer " + \
            "-i " + \
            str(infile) + \
            " " + \
            "-o " + \
            str(outfile) + \
            " " + \
            "--tglobal=" + \
            str(shift)
        run_cmd = main_cmd
        subprocess.run(run_cmd, shell=True)

    temp_shift(infile=str(src_prefix)+"/odd.nii.gz",
               outfile=str(src_prefix)+"/odd_earlier.nii.gz",
               shift=0.25)
    temp_shift(infile=str(src_prefix)+"/odd.nii.gz",
               outfile=str(src_prefix)+"/odd_later.nii.gz",
               shift=-0.25)
    temp_shift(infile=str(src_prefix)+"/even.nii.gz",
               outfile=str(src_prefix)+"/even_earlier.nii.gz",
               shift=0.25)
    temp_shift(infile=str(src_prefix)+"/even.nii.gz",
               outfile=str(src_prefix)+"/even_later.nii.gz",
               shift=0.75)

    main_cmd = "fslmaths " + \
        str(src_prefix)+"/odd_later.nii.gz" + \
        " " + \
        "-sub " + \
        str(src_prefix)+"/even_earlier.nii.gz" + \
        " " + \
        str(src_prefix)+"/even_to_odd.nii.gz"
    run_cmd = main_cmd
    subprocess.run(run_cmd, shell=True)

    main_cmd = "fslmaths " + \
        str(src_prefix)+"/odd_earlier.nii.gz" + \
        " " + \
        "-sub " + \
        str(src_prefix)+"/even_later.nii.gz" + \
        " " + \
        str(src_prefix)+"/odd_to_even.nii.gz"
    run_cmd = main_cmd
    subprocess.run(run_cmd, shell=True)

    surr_diff_calc_path = Path(str(src_prefix)+"/surr_diff")
    make_dir(src=surr_diff_calc_path)

    even_sub_path = Path(str(surr_diff_calc_path)+"/e_")
    split4d(src=Path(str(src_prefix)+"/even_to_odd.nii.gz"),
            out=even_sub_path)
    even_sub_filelist = list(surr_diff_calc_path.glob("e_*.nii.gz"))
    even_sub_filelist.sort()

    odd_sub_path = Path(str(surr_diff_calc_path)+"/o_")
    split4d(src=Path(str(src_prefix)+"/odd_to_even.nii.gz"),
            out=odd_sub_path)
    odd_sub_filelist = list(surr_diff_calc_path.glob("o_*.nii.gz"))
    odd_sub_filelist.sort()

    subtracted_fulllist = even_sub_filelist + odd_sub_filelist
    subtracted_fulllist[::2] = even_sub_filelist
    subtracted_fulllist[1::2] = odd_sub_filelist

    subtracted_fulllist_merge = list_to_string(in_list=subtracted_fulllist)
    merge3d(src=subtracted_fulllist_merge,
            out=str(src_prefix)+"_moCorr_surrDiff.nii.gz",
            tr=2.861)
    calc_Tmean(src=Path(str(src_prefix)+"_moCorr_surrDiff.nii.gz"))

    if remean is True:
        mean_out = calc_Tmean(src)
        main_cmd = "fslmaths " + \
            str(src_prefix)+"_moCorr_surrDiff.nii.gz" + \
            " " + \
            "-add " + \
            str(src_prefix)+"_Tmean.nii.gz" + \
            " " + \
            str(src_prefix)+"_moCorr_surrDiff_plusMean.nii.gz"
        run_cmd = main_cmd
        subprocess.run(run_cmd, shell=True)

    remove_dir(surr_diff_calc_path)

    if surr_avg is True:
        main_cmd = "fslmaths " + \
            str(src_prefix)+"/odd_later.nii.gz" + \
            " " + \
            "-add " + \
            str(src_prefix)+"/even_earlier.nii.gz" + \
            " " + \
            "-div 2 " + \
            str(src_prefix)+"/even_to_odd.nii.gz"
        run_cmd = main_cmd
        subprocess.run(run_cmd, shell=True)

        main_cmd = "fslmaths " + \
            str(src_prefix)+"/odd_earlier.nii.gz" + \
            " " + \
            "-add " + \
            str(src_prefix)+"/even_later.nii.gz" + \
            " " + \
            "-div 2 " + \
            str(src_prefix)+"/odd_to_even.nii.gz"
        run_cmd = main_cmd
        subprocess.run(run_cmd, shell=True)

        surr_avg_calc_path = Path(str(src_prefix)+"/surr_avg")
        make_dir(src=surr_avg_calc_path)

        even_sub_path = Path(str(surr_avg_calc_path)+"/e_")
        split4d(src=Path(str(src_prefix)+"/even_to_odd.nii.gz"),
                out=even_sub_path)
        even_sub_filelist = list(surr_avg_calc_path.glob("e_*.nii.gz"))
        even_sub_filelist.sort()

        odd_sub_path = Path(str(surr_avg_calc_path)+"/o_")
        split4d(src=Path(str(src_prefix)+"/odd_to_even.nii.gz"),
                out=odd_sub_path)
        odd_sub_filelist = list(surr_avg_calc_path.glob("o_*.nii.gz"))
        odd_sub_filelist.sort()

        averaged_fulllist = even_sub_filelist + odd_sub_filelist
        averaged_fulllist[::2] = even_sub_filelist
        averaged_fulllist[1::2] = odd_sub_filelist

        averaged_fulllist_merge = list_to_string(in_list=averaged_fulllist)
        merge3d(src=averaged_fulllist_merge,
                out=str(src_prefix)+"_moCorr_surrAvg.nii.gz",
                tr=2.861)
        calc_Tmean(src=Path(str(src_prefix)+"_moCorr_surrAvg.nii.gz"))

        remove_dir(surr_avg_calc_path)

    remove_files(src=src_prefix,
                 pattern="vol*.nii.gz")
    remove_files(src=src_prefix,
                 pattern="odd*.nii.gz")
    remove_files(src=src_prefix,
                 pattern="even*.nii.gz")

def despike_fmap(fmap,fmap_mask):
    """[summary]

    Args:
        fmap ([type]): [description]
        fmap_mask ([type]): [description]
    """
    if type(fmap) is str:
        fmap = Path(fmap)
    if type(fmap_mask) is str:
        fmap_mask = Path(fmap_mask)

    fmap_prefix = get_prefix(fmap)
    fmap_mask_prefix = get_prefix(fmap_mask)
    
    main_cmd = "fslmaths " + \
        str(fmap_mask) + \
        " " + \
        "-ero " + \
        str(fmap_mask_prefix)+"_ero.nii.gz"
    run_cmd = main_cmd
    subprocess.run(run_cmd, shell=True)
    
    main_cmd = "fugue " + \
        "--loadfmap=" + \
        str(fmap) + \
        " " + \
        "--savefmap=" + \
        str(fmap_prefix)+"_filt.nii.gz" + \
        " " + \
        "--mask=" + \
        str(fmap_mask_prefix)+"_ero.nii.gz" + \
        " "
    post_cmd = "--despike " + \
        "--despikethreshold=2.1"
    run_cmd = main_cmd+post_cmd
    subprocess.run(run_cmd, shell=True)

    main_cmd = "fslmaths " + \
        str(fmap) + \
        " " + \
        "-sub " + \
        str(fmap_prefix)+"_filt.nii.gz" + \
        " " + \
        "-mas " + \
        str(fmap_mask) + \
        " " + \
        "-add " + \
        str(fmap_prefix)+"_filt.nii.gz" + \
        " " + \
        str(fmap_prefix)+"_despiked.nii.gz"
    run_cmd = main_cmd
    subprocess.run(run_cmd, shell=True)
    
    file_rm=Path(str(fmap_prefix)+"_filt.nii.gz")
    file_rm.unlink()
    file_rm=Path(str(fmap_mask_prefix)+"_ero.nii.gz")
    file_rm.unlink()

def flirt_transform(ref_file_path, trg_file_path, mat_file_path, verbose=False, sim=False):
    """Apply estimated transformations using flirt

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

    main_cmd = "flirt -interp sinc -applyxfm " + \
            "-init " + mat_file_path.as_posix() + " " + \
            "-ref " + ref_file_path.as_posix() + " " + \
            "-in " + trg_file_path.as_posix() + " " + \
            "-out " + out_file_prefix.as_posix() + "_flirtWarped.nii.gz"

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

def flirt_bbregister(rage_file_path,wm_seg_path,tse_file_path,mat_file_path,verbose=False):
    """Register the TSE to hires MP2RAGE using an init matrix & BBR,
       invert and apply to hires MP2RAGE to transform to TSE space

    Args:
        rage_file_path ([type]): Path to hires MP2RAGE file
        wm_seg_path ([type]): Path to WM Seg from hires MP2RAGE file
        tse_file_path ([type]): Path to T2w TSE file
        mat_file_path ([type]): Path to init mat file
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
    
    main_cmd_0 = "flirt -interp sinc -nosearch -noresample -noresampblur " + \
        "-schedule ${FSLDIR}/etc/flirtsch/bbr.sch -cost bbr -minsampling 0.3 " + \
        "-init " + mat_file_path.as_posix() + " " + \
        "-ref " + rage_file_path.as_posix() + " " + \
        "-wmseg " + wm_seg_path.as_posix() + " " + \
        "-in " + tse_file_path.as_posix() + " " + \
        "-omat " + out_file_prefix.as_posix() + "_fsl_bbReg.mat"  + " " + \
        "-out " + out_file_prefix.as_posix() + "_fsl_bbReg.nii.gz"

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
    
    main_cmd_1 = "convert_xfm " + \
    "-omat " + out_file_prefix.as_posix() + "_fsl_bbReg_inverse.mat" + \
    " " + \
    "-inverse " + out_file_prefix.as_posix() + "_fsl_bbReg.mat"

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
    
    flirt_transform(ref_file_path=tse_file_path,
                    trg_file_path=rage_file_path,
                    mat_file_path=Path(out_file_prefix.as_posix() + "_fsl_bbReg_inverse.mat"))
    
    