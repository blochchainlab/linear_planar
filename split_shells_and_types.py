import nibabel as nib
import numpy as np
import sys


def main(dwifile, bvals, gtablefile, bmax, eddy_configfile, eddy_indicesfile, dwiout, bvecout, bvalout, eddy_indices_out, position_idx_out):

	# convert cmdline strings
	bvals = [int(i) for i in bvals.split(',')]
	bmax = int(bmax)

	# load data
	dwiimg = nib.load(dwifile)
	dwiall = dwiimg.get_data()

	gtable = np.genfromtxt(gtablefile)
	eddy_indices = np.genfromtxt(eddy_indicesfile)

	# guess bval and type from gradient table
	gnorm = np.linalg.norm(gtable[:, :3], axis=1)
	# assumes g1 and g2 have same bvalue
	databval = bmax*gnorm**2

	# linear-planar detection (assumes only perfect linear and planar exist (eucl. dist. between g1 and g2 is 0 OR sqrt(2)~1.4))
	isPlanar = np.linalg.norm((gtable[:, :3]-gtable[:, 3:6])/gnorm[:,None], axis=1) > 1
	# detect which shell in bvals fits best dataval
	shell_idx = np.argmin(np.abs(databval - np.array(bvals)[:,None]), axis=0)

	# split it
	for idx, bval in enumerate(bvals):
		# checks for linear, if any
		b_mask = np.logical_and(shell_idx==idx, isPlanar==0)
		if np.sum(b_mask) > 0:
			print('{} linear b = {} volumes'.format(np.sum(b_mask), bval))
			position_idx = np.where(b_mask)[0]
			np.savetxt(position_idx_out + '__b{}_lin.txt'.format(bval), position_idx, fmt='%i')
			np.savetxt(eddy_indices_out + '__b{}_lin.txt'.format(bval), eddy_indices[position_idx], fmt='%i')
			np.savetxt(bvalout + '__b{}_lin.txt'.format(bval), databval[position_idx])
			dirs = gtable[position_idx, :3]
			normdirs = dirs / np.linalg.norm(dirs, axis=1)[:,None]
			np.savetxt(bvecout + '__b{}_lin.txt'.format(bval), normdirs)
			newimg = nib.nifti1.Nifti1Image(dwiall[...,position_idx], dwiimg.affine, dwiimg.header)
			nib.save(newimg, dwiout + '__b{}_lin.nii'.format(bval))
			newimg2 = nib.nifti1.Nifti1Image(dwiall[...,position_idx].mean(3), dwiimg.affine, dwiimg.header)
			nib.save(newimg2, dwiout + '__b{}_lin_SM.nii'.format(bval))

		# checks for planar, if any
		b_mask = np.logical_and(shell_idx==idx, isPlanar==1)
		if np.sum(b_mask) > 0:
			print('{} planar b = {} volumes'.format(np.sum(b_mask), bval))
			position_idx = np.where(b_mask)[0]
			np.savetxt(position_idx_out + '__b{}_pla.txt'.format(bval), position_idx, fmt='%i')
			np.savetxt(eddy_indices_out + '__b{}_pla.txt'.format(bval), eddy_indices[position_idx], fmt='%i')
			np.savetxt(bvalout + '__b{}_pla.txt'.format(bval), databval[position_idx])
			dirs = gtable[position_idx, :3] + gtable[position_idx, 3:6]
			normdirs = dirs / np.linalg.norm(dirs, axis=1)[:,None]
			np.savetxt(bvecout + '__b{}_pla.txt'.format(bval), normdirs)
			newimg = nib.nifti1.Nifti1Image(dwiall[...,position_idx], dwiimg.affine, dwiimg.header)
			nib.save(newimg, dwiout + '__b{}_pla.nii'.format(bval))
			newimg2 = nib.nifti1.Nifti1Image(dwiall[...,position_idx].mean(3), dwiimg.affine, dwiimg.header)
			nib.save(newimg2, dwiout + '__b{}_pla_SM.nii'.format(bval))

if __name__ == "__main__":
	# for idx,i in enumerate(sys.argv):
	# 	print(idx, i)
	main(*sys.argv[1:])

