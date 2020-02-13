import numpy as np
import nibabel as nib
import scipy.optimize as opt
# import pylab as pl

from smt_lin_pla import *

from multiprocess import Pool

import sys
import os

from time import time



## for multiprocess
N_pool = 6







## This part would be kind-of how to create b and bd from G1G2G3.txt instead of hardcoding
# gtable = np.loadtxt('/home/raid2/paquette/work/genDirsPlanar/scheme/test_45_45_x3.txt')
# ten = table2ten(gtable, bval=1.5)
# # convert s/mm^2 to ms/um^2
# # ten = btensor / 1000.
# buni, bval = discritizeBvals(ten)
# scheme, scheme_norm, dir_type = bd_getDirs(ten)








## load data
masterpath = '/some/path/'

# in this case, we had a b50 instead of a b0, so I create a fake_B0 by extrapolating from a ensor fit on the b50+b1000, see preprocessing.sh
b0_path = masterpath + 'PREPRO/computed_S0.nii'

# this structure assumes that each bvalue and b-shape are save separately in already spherical averaged form
data_path = [masterpath + 'PREPRO/v1_eddy_b1000_lin_SM.nii.gz',
			 masterpath + 'PREPRO/v1_eddy_b1000_pla_SM.nii.gz',
			 masterpath + 'PREPRO/v1_eddy_b2000_lin_SM.nii.gz',
			 masterpath + 'PREPRO/v1_eddy_b2000_pla_SM.nii.gz',
			 masterpath + 'PREPRO/v1_eddy_b4000_lin_SM.nii.gz',
			 masterpath + 'PREPRO/v1_eddy_b4000_pla_SM.nii.gz']

mask_path = masterpath + 'PREPRO/bet_init_b50_dil.nii.gz'

savepath = masterpath + 'FIT/fit1/'


## hardcoded b-tensor size and shape corresponding to the load order of data_path (in um^2/ms)
b = np.array([1., 1., 2., 2., 4., 4.])
bd = np.array([1., -.5, 1., -.5, 1., -.5])


if not os.path.isdir(savepath):
	print('saving directory doesn\'t exist')
	# return None


b0_ima = nib.load(b0_path)
b0_data = b0_ima.get_data()
data = np.concatenate([nib.load(path).get_data()[...,None] for path in data_path], axis=3)
mask = nib.load(mask_path).get_data().astype(bool)
print('data loaded')




## data b0 normalization, clip outlier to [0, 1]
data_norm = np.clip(data, 0, np.inf) / np.clip(b0_data, 0, np.inf)[...,None]
data_norm[np.isnan(data_norm)] = 0
# clip outlier
data_SM = np.clip(data_norm, 0, 1)
print('data normalized')





## defining objective function from signal function of smt_lin_pla (model2)
## note: if using model 1, f_csf is no longuer a parameter

# single-voxel per-shell error function
def func_vox(p, b, bd, datavox):
	# unpack param
	Dpar_in, Dpar_ex, Dperp_ex, f_in, f_ex, f_csf = p
	# return model(param)-data for model2
	return sm_signal_model_2(b, bd, Dpar_in, Dpar_ex, Dperp_ex, f_in, f_ex, f_csf) - datavox

# single-voxel squared_error 
def func_vox_ls(p, b, bd, datavox):
	# compute single-voxel per-shell error function
	f = func_vox(p, b, bd, datavox)
	# return  sum_{shell} (model(param)-data)^2
	return np.sum(f**2, axis=0)



#################################################
## CONSTRAINT set 1, to use when NOT USING the ratio non-linear constraint
## note: if using model 1, 
## the equality constraint 1 <= f_in + f_ex + f_csf <= 1 
## would become the inequality constraint 0 <= f_in + f_ex <= 1 

# # 0.1 <= D_in <= 3
# # 0.1 <= D_ex_par <= 3
# # 0.1 <= D_ex_perp <= 1.5
# # 0 <= f_in <= 1
# # 0 <= f_ex <= 1
# # 0 <= f_csf <= 1
# # 0 <= D_ex_par - D_ex_perp <= 3
# # 1 <= f_in + f_ex + f_csf <= 1

# A = np.array([[1, 0, 0, 0, 0, 0],
# 			  [0, 1, 0, 0, 0, 0],
# 			  [0, 0, 1, 0, 0, 0],
# 			  [0, 0, 0, 1, 0, 0],
# 			  [0, 0, 0, 0, 1, 0],
# 			  [0, 0, 0, 0, 0, 1],
# 			  [0, 1, -1, 0, 0, 0],
# 			  [0, 0, 0, 1, 1, 1]])

# lb = np.array([0.1,0.1,0.1,0,0,0,0,1])
# ub = np.array([3,3,1.5,1,1,1,3,1])
#################################################


#################################################
## CONSTRAINT set 2, to use when USING the ratio non-linear constraint
## note: if using model 1, 
## the equality constraint 1 <= f_in + f_ex + f_csf <= 1 
## would become the inequality constraint 0 <= f_in + f_ex <= 1

# because 0 <= D_ex_par / D_ex_perp <= 6 in the nonlin constraint, we can simplify
# 0.1 <= D_in <= 3
# 0.1 <= D_ex_par <= 3
# 0.1 <= D_ex_perp <= 1.5
# 0 <= f_in <= 1
# 0 <= f_ex <= 1
# 0 <= f_csf <= 1
# 1 <= f_in + f_ex + f_csf <= 1

A = np.array([[1, 0, 0, 0, 0, 0],
			  [0, 1, 0, 0, 0, 0],
			  [0, 0, 1, 0, 0, 0],
			  [0, 0, 0, 1, 0, 0],
			  [0, 0, 0, 0, 1, 0],
			  [0, 0, 0, 0, 0, 1],
			  [0, 0, 0, 1, 1, 1]])

lb = np.array([0.1,0.1,0.1,0,0,0,1])
ub = np.array([3,3,1.5,1,1,1,1])
#################################################



## linear constraint object
cons = opt.LinearConstraint(A, lb, ub)



#################################################
## non-linear ratio constraint set

## constraint the Dpar / Dperp ratio of the extra compartement
## equivalent to a constraint on the FA of the compartement

## unclear how necessary it is
## pretty time-costly


def ratio(x):
	# Dpar / Dperp
	return x[1]/x[2]

# derivative of the constraint for all parameters
def ratio_jac(x):
	# d/dDin = 0
	# d/dDpar = 1/Dperp
	d1 = 1 / x[2]
	# d/dPerp = -Dpar / Dperp^2
	d2 = -x[1] / x[2]**2
	# d/dfin = 0
	# d/dfex = 0
	# d/dfw = 0
	return np.array([0,d1,d2,0,0,0])

## some FA to ratio conversion, computed elsewhere
# # FA ~ 0.2425
# ratio_min = 1.5

# FA ~ 0.0
ratio_min = 1.

# FA ~ 0.8111
ratio_max = 6.0
#################################################


## nonlinear ratio constraint object
fa_cons = opt.NonlinearConstraint(ratio, np.array([ratio_min]), np.array([ratio_max]), jac=ratio_jac, hess=opt.BFGS())






## uFA for this type of models
def microFA(d1, d2, d3, fin, fex, fcsf):
	# if (fin < 0) or (fex < 0) or (fin + fex > 1):
	# 	return np.nan

	microAx = fin*d1 + fex*d2 + fcsf*3
	microRad = fex*d3 + fcsf*3
	microMean = (microAx + 2*microRad) / 3.
	microFA = np.sqrt((3/2.)*((microAx - microMean)**2 + 2*(microRad - microMean)**2) / (microAx**2 + 2*microRad**2))
	return microFA


## compartement FA
def FA(Dpar, Dperp):
	meanD = (Dpar+2*Dperp)/3.
	return np.sqrt(3/2.) * np.sqrt((Dpar-meanD)**2 + 2*(Dperp-meanD)**2) / np.sqrt(Dpar**2 + 2*Dperp**2)







## Start point for the optimization
## practice seem to indicate this is not super sensitive
## by principle, i'm happier when I initialize my optimization away from my constraint boundaries
x0 = np.array([2.0, 1.3, 0.7, 0.65, 0.3, 0.05])




## storage
function_residual = np.zeros(data_SM.shape[:3])
exit_code = np.zeros(data_SM.shape[:3])
fitted_param = np.zeros(data_SM.shape[:3]+(6,)) ## this would be (5,) for model1



# voxelwise minimization with linear and ratio constraint
def run_optim_fixed_x0(datavox):
	res = opt.minimize(func_vox_ls, x0=x0, args=(b, bd, datavox), jac='3-point', hess=opt.BFGS(), method='trust-constr', constraints=[cons, fa_cons])
	return res



## loop over voxels in the mask with multiprocess.Pool.map
print('setting up pool ({})'.format(N_pool))

datavoxs = data_SM[mask]

pool = Pool(N_pool)

t1 = time()
res = pool.map(run_optim_fixed_x0, datavoxs)
t2 = time()

print('time = {} s\n'.format(t2 - t1))


## store values
function_residual[mask] = np.array([r.fun for r in res])
exit_code[mask] = np.array([r.status for r in res])
fitted_param[mask] = np.array([r.x for r in res])


## save
endtag = '1'
nib.save(nib.nifti1.Nifti1Image(fitted_param, b0_ima.affine, b0_ima.header), savepath + 'fit{}.nii.gz'.format(endtag))
nib.save(nib.nifti1.Nifti1Image(exit_code, b0_ima.affine, b0_ima.header), savepath + 'exit{}.nii.gz'.format(endtag))
nib.save(nib.nifti1.Nifti1Image(function_residual, b0_ima.affine, b0_ima.header), savepath + 'res{}.nii.gz'.format(endtag))


## compute uFA voxelwise and save
uFA = np.zeros(data_SM.shape[:3])
for vox in np.ndindex(data_SM.shape[:3]):
	if mask[vox]:
		uFA[vox] = microFA(*fitted_param[vox])

nib.save(nib.nifti1.Nifti1Image(uFA, b0_ima.affine, b0_ima.header), savepath + 'uFA{}.nii.gz'.format(endtag))



## compute FA of extra compartement (sanity check /  necessity check of nonlin constraint) voxelwise and save

FAextra = np.zeros(data_SM.shape[:3])
for vox in np.ndindex(data_SM.shape[:3]):
	if mask[vox]:
		FAextra[vox] = FA(fitted_param[vox][1],fitted_param[vox][2])


nib.save(nib.nifti1.Nifti1Image(FAextra, b0_ima.affine, b0_ima.header), savepath + 'FA_extra{}.nii.gz'.format(endtag))


print('Done')
pool.close()




