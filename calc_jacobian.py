#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import nibabel as nib
import sys
import argparse
import os

# Hard-coded max jacobian determinant for siemens, taken from grad unwarp file of Human Connectom Project
siemens_max_det = 10.

DESCRIPTION =   'Calculation of the Jacobi Determinant from a given warp field. Cornelius Eichner 2018'

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('--in', dest='data', action='store', type=str,
                            help='Path of the input warp field (nifti format)')

    p.add_argument('--out', dest='out', action='store', type=str,
                            help='Path of the output')

    return p


def evaluate_jacobian(F, pixdim):

    # Differentials in each dimension
    d0, d1, d2 = pixdim[0], pixdim[1], pixdim[2]

    if d0 == 0 or d1 == 0 or d2 == 0:
        raise ValueError('weirdness found in Jacobian calculation')

    dFxdx, dFxdy, dFxdz = np.gradient(F[...,0], d0, d1, d2)
    dFydx, dFydy, dFydz = np.gradient(F[...,1], d0, d1, d2)
    dFzdx, dFzdy, dFzdz = np.gradient(F[...,2], d0, d1, d2)

    jacdet = (dFxdx) * (dFydy) * (dFzdz) \
           - (dFxdx) * dFydz * dFzdy \
           - dFxdy * dFydx * (dFzdz) \
           + dFxdy * dFydz * dFzdx \
           + dFxdz * dFydx * dFzdy \
           - dFxdz * (dFydy) * dFzdx

    jacdet = np.abs(jacdet)
    jacdet[np.where(jacdet > siemens_max_det)] = siemens_max_det

    return jacdet



def main():
    # Load parser to read data from command line input
    parser = buildArgsParser()
    args = parser.parse_args()

    # Change path to location of input data
    # PATH = os.path.dirname(os.path.realpath(args.data)) + '/'
    # os.chdir(PATH)

    DATA = os.path.realpath(args.data)
    OUT = args.out

    # Read data and extract necessary information
    field = nib.load(DATA)
    field_data = field.get_data()
    hd_pixdim = field.header.get('pixdim')
    pixdim = hd_pixdim[1:4]

    jacdet = evaluate_jacobian(field_data, pixdim)

    nib.save(nib.Nifti1Image(jacdet, field.affine), OUT)


if __name__ == '__main__':
    main()
