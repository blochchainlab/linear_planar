#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import nibabel as nib
import numpy as np

DESCRIPTION =   'Joint correction of eddy currents and gradient non linearity correction using FSL tools. Cornelius Eichner 2018'

PATH_GRAD_UNWARP = '/home/raid2/paquette/tools/gradunwarp/gradunwarp/core/gradient_unwarp.py'
# PATH_EDDY = '/home/raid2/paquette/tools/eddy_openmp'
PATH_EDDY = '/home/raid2/paquette/tools/eddy_cuda8.0'
PATH_CALC_JACOBIAN = '/home/raid2/paquette/tools/onestep_eddy_nlgc/calc_jacobian.py'


def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('--in', dest='data', action='store', type=str,
                            help='Path of the input volume (nifti format)')

    p.add_argument('--bvec', dest='bvec', action='store', type=str,
                            help='Path of the bvec file (default \'bvec\')')

    p.add_argument('--bval', dest='bval', action='store', type=str,
                           help='Path of the bval file (default \'bval\')')

    p.add_argument('--mask', dest='mask', action='store', type=str,
                            help='Path of the brain mask')

    p.add_argument('--acqp', dest='acqp', action='store', type=str,
                        help='Path of eddy acquisition parameter file (default \'acq_param.txt\')')
    
    p.add_argument('--index', dest='index', action='store', type=str,
                            help='Path of the eddy index file (default \'index.txt\')')

    p.add_argument('--topup', dest='topup', action='store', type=str,
                           help='Base name of topup file structure (default \'topup/topup\')')

    p.add_argument('--out', dest='out', action='store', type=str,
                            help='Path of the output volume')

    p.add_argument('--pmat', dest='pmat', action='store', type=str,
                            help='Path of the affine transform to apply at the end of registrationb (e.g., T1 anatomy)')

    p.add_argument('--mp', dest='openmp', action='store', type=int, default='1', 
                            help='Optional: OpenMP parallelization factor (default 8)')

    p.add_argument('--jump', dest='jump', action='store', type=str,
                           help='skip part of the processing')

    return p


def main():
    ######################
    # Read Input Arguments

    parser = buildArgsParser()
    args = parser.parse_args()

    if args.bvec is None:
        BVEC = 'bvec'
    else:
        BVEC = os.path.realpath(args.bvec)

    if args.bval is None:
        BVAL = 'bval'
    else:
        BVAL = os.path.realpath(args.bval)

    if args.acqp is None:
        ACQP = 'acq_param.txt'
    else:
        ACQP = os.path.realpath(args.acqp)

    if args.index is None:
        INDEX = 'index.txt'
    else:
        INDEX = os.path.realpath(args.index)

    if args.topup is None:
        TOPUP = 'topup/topup'
    else:
        TOPUP = os.path.realpath(args.topup)

    if args.openmp is None:
        OPENMP_THREADS = 8
    else:
        OPENMP_THREADS = args.openmp

    if args.pmat is None:
        PMAT = None
    else:
        PMAT = args.pmat

    if args.jump is None:
        jump = 0
    else:
        jump = int(args.jump)

    ######################
    # Change Pathways to operate in data folder

    PATH = os.path.dirname(os.path.realpath(args.data)) + '/'
    os.chdir(PATH)

    DATA = os.path.realpath(args.data)
    MASK = os.path.realpath(args.mask)
    
    OUT = args.out


    ######################
    # Calculate Gradient Non Linearities
    if jump < 1:
        os.system("echo Correction of gradient non linearities")

        # Extract first data volume for correction
        cmd = 'fslroi ' + DATA + ' ' + PATH + 'single.nii.gz 0 1'
        os.system(cmd)

        # Perform gradient non linearity correction
        cmd = PATH_GRAD_UNWARP + ' '+ PATH + 'single.nii.gz ' + PATH + 'nlgc.nii.gz siemens -g /home/raid2/paquette/tools/gradunwarp/coeffConnectom.grad -n'
        os.system(cmd)
        
        # Rename resulting files
        cmd = 'mv fullWarp_abs.nii.gz ' + PATH + 'nlgc_warp.nii.gz'
        os.system(cmd)
    else:
        print('jumping over nlgc warp estimation')

    """
    Resulting output files:
    single.nii.gz           first volume of data extracted, as gradient unwarper cannot handle 4D data
    nlgc.nii.gz             first volume, corrected for gradient non linearities, will not be further employed
    nlgc_warp.nii.gz        warp field for gradient correction
    """


    ######################
    # Calculate eddy displacement fields
    if jump < 2:
        os.system("echo Running FSL Eddy ")

        # Set number of parallel processing threads
        cmd = 'set OMP_NUM_THREADS = ' + str(OPENMP_THREADS)
        os.system(cmd)

        # Run eddy from home installation
        cmd =\
        PATH_EDDY + ' \
            --imain='+ DATA + ' \
            --mask=' + MASK + ' \
            --index=' + INDEX + ' \
            --acqp=' + ACQP + ' \
            --bvecs=' + BVEC + ' \
            --bvals=' + BVAL + ' \
            --topup=' + TOPUP + ' \
            --out=' + PATH + 'eddy_tmp \
            --dfields= \
            --repol \
            --data_is_shelled \
            --cnr_maps'

        os.system(cmd)
        os.system("echo Eddy Complete")


        # Move all displacement field files in dfield subdirectory
        os.system("echo Moving displacement field files into common folder")
        os.system('mkdir -p dfield')
        os.system('mv *displacement_fields* dfield/')
    else:
        print('jumping over eddy')
    DFIELD_FILES = sorted(os.listdir('dfield/'))

    """
    Resulting output files:
    eddy_tmp.nii.gz         Output of FSL eddy, data will not be further employed since warps were already applied
    eddy_tmp.eddy_outlier_free_data.nii.gz
                            Outlier free data, from FSL eddy without potential signal dropouts - file will be employed for volume warping
    eddy_tmp.eddy_rotated_bvecs
                            Rotated bvecs file from eddy correction
    dfield/                 Folder containing all final deformation fields from eddy
    """


    ######################
    # Combine both warp fields 
    if jump < 3:
        os.system("echo Combining Warp Fields")

        os.system('mkdir -p comb_warp')

        # Combination of warpfields by concatenation
        if PMAT is None:
            for i in DFIELD_FILES:
                cmd = 'convertwarp \
                            -o comb_warp/nlcg.' + i + '\
                            -r nlgc_warp.nii.gz \
                            --warp1=dfield/' + i + ' \
                            --warp2=nlgc_warp.nii.gz '
                os.system(cmd)
        else:
            for i in DFIELD_FILES:
                cmd = 'convertwarp \
                            -o comb_warp/nlcg.' + i + ' \
                            -r nlgc_warp.nii.gz \
                            --postmat=' + PMAT + ' \
                            --warp1=dfield/' + i + ' \
                            --warp2=nlgc_warp.nii.gz '
                os.system(cmd)
    else:
        print('jumping over combine warp')

    COMB_WARP_FILES = sorted(os.listdir('comb_warp/'))

    """
    Resulting output files:
    comb_warp/              Folder containing all combined warp files    
    """


    ######################
    # Apply warp for using combined field
    os.system("echo Application of combined warp fields")

    os.system('mkdir -p split_data')
    os.system('mkdir -p warped_data')

    # cmd = 'fslsplit ' + PATH + '*eddy_outlier_free*.nii.gz split_data/data -t'
    cmd = 'fslsplit ' + PATH + 'eddy_tmp.eddy_outlier_free_data.nii.gz split_data/data -t'
    os.system(cmd)

    SPLIT_DATA_FILES = sorted(os.listdir('split_data/'))

    if jump < 4:
        for i in xrange(len(os.listdir('split_data/'))):
            cmd = 'applywarp \
                    -i split_data/' + SPLIT_DATA_FILES[i] + ' \
                    -r split_data/' + SPLIT_DATA_FILES[0] + ' \
                    -o warped_data/' + SPLIT_DATA_FILES[i] + ' \
                    -w comb_warp/' + COMB_WARP_FILES[i] + ' \
                    --interp=spline \
                    --datatype=float'
            os.system(cmd)
    else:
        print('jumping over applywarp')


    ######################
    # Correct for warp signal bias using Jacobian determinante
    os.system("echo Correction of signal intensities with Jacobian Determinant")

    os.system('mkdir -p jac_det_data')
    os.system('mkdir -p corr_data')

    for i in xrange(len(os.listdir('comb_warp/'))):
        cmd = 'python ' + PATH_CALC_JACOBIAN + ' --in comb_warp/' + COMB_WARP_FILES[i] + ' --out jac_det_data/jac.' + COMB_WARP_FILES[i]
        os.system(cmd)
        cmd = 'fslmaths warped_data/' + SPLIT_DATA_FILES[i] + ' -mul jac_det_data/jac.' + COMB_WARP_FILES[i] + ' corr_data/' + SPLIT_DATA_FILES[i]
        # cmd = 'fslmaths warped_data/' + SPLIT_DATA_FILES[i] + ' -div jac_det_data/jac.' + COMB_WARP_FILES[i] + ' corr_data/' + SPLIT_DATA_FILES[i]
        os.system(cmd)


    ######################
    # Merge warped volumes
    os.system("echo Merging corrected volumes")

    cmd = 'fslmerge -t ' + OUT + ' corr_data/*'
    os.system(cmd)

    """
    Resulting output files:
    OUT                     Final warped output file
    split_data/             Folder containing all volumes as single files, necessary for volume warping    
    """


if __name__ == '__main__':
    main()
