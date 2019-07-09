import nibabel as nib
import numpy as np
import sys

if __name__ == "__main__":
	print('\nparam order::\ninput_vol.nii(.gz) output_vol_mean.nii(.gz) output_vol_std.nii(.gz)\n')

	# loading data
	img = nib.load(sys.argv[1])
	data = img.get_data()

	newdata1 = data.mean(axis=3)
	newdata2 = data.std(axis=3)

	newimg1 = nib.nifti1.Nifti1Image(newdata1, img.affine, img.header)
	newimg2 = nib.nifti1.Nifti1Image(newdata2, img.affine, img.header)

	nib.save(newimg1, sys.argv[2])
	nib.save(newimg2, sys.argv[3])


