# basic hardcoded Dicom to Nifti script using mrtrix-suite


FOLDER=/path/to/data/


# non-diffusion data with "topup" export flag
# the stride flag is extremely important
mrconvert $FOLDER'dicoms/' $FOLDER'data.nii' -stride 1,2,3 -datatype float64 -nthreads 6 -export_pe_table $FOLDER'PE.txt'


# diffusion data with "topup+eddy" export flag and bvec/bval export flag
# the stride flag is extremely important
mrconvert $FOLDER'dicoms/' $FOLDER'data.nii' -stride 1,2,3,4 -datatype float64 -nthreads 6 -export_pe_eddy $FOLDER'PE.txt' $FOLDER'eddy_indices.txt' --export_grad_fsl $FOLDER'bvec.txt' $FOLDER'bval.txt'

