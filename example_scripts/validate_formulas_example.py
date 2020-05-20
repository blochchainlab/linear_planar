import numpy as np
from gs import findLinIndepRandomRot, gramSchmidt3
from dipy.core.sphere import HemiSphere, disperse_charges
# from sphHist import plotScatter3
from smt_lin_pla import *


## we want to generate a set of well distribution point on the sphere
## and then generate a linear and a planar B-tensor for each point
## the linear one points in the generated direction
## the planar one have their plane normal pointing toward it

## This is a bit old, it uses some naming convention that are different then the fitting code


## Generate semi-well distributed point on the sphere with golden spiral method 
N = 1000
golden_angle = np.pi * (3 - np.sqrt(5))
theta = golden_angle * np.arange(N)
z = np.linspace(1 - 1.0 / N, 1.0 / N - 1, N)
radius = np.sqrt(1 - z * z)
 
points = np.zeros((N, 3))
points[:,0] = radius * np.cos(theta)
points[:,1] = radius * np.sin(theta)
points[:,2] = z


## convert to dipy hemisphere
hemi = HemiSphere(xyz=points)
res, _ = disperse_charges(hemi, 100)
dirs = res.vertices


## get necessary vector to generate the plane
planar_dirs = np.zeros((dirs.shape[0], 2, 3))
for idx in range(dirs.shape[0]):
	# current directions
	v = dirs[idx]
	# generate initial non colinear vectors, v1 is v
	v1, v2, v3 = findLinIndepRandomRot(v)
	# generate 3 orthonormal vectors, u1 is v normalized
	u1, u2, u3 = gramSchmidt3(v1, v2, v3)
	planar_dirs[idx, 0] = u2
	planar_dirs[idx, 1] = u3


## G1G2G3.txt-like structure
## g3 = 0 because we are not doing spherical B-tensor
g1 = np.concatenate((dirs, planar_dirs[:,0]))
g2 = np.concatenate((dirs, planar_dirs[:,1]))
g3 = np.zeros_like(g1)
grad = np.concatenate((g1,g2,g3), axis=1)


## 3-gradient form to B-tensor
def table2ten(gtable, bval = 1000.):
	tensors_1 = np.zeros((gtable.shape[0], 3, 3))
	tensors_2 = np.zeros((gtable.shape[0], 3, 3))
	tensors_3 = np.zeros((gtable.shape[0], 3, 3))

	for ix in range(3):
		for iy in range(3):
			tensors_1[:, ix, iy] = gtable[:, ix]*gtable[:, iy]
			tensors_2[:, ix, iy] = gtable[:, ix+3]*gtable[:, iy+3]
			tensors_3[:, ix, iy] = gtable[:, ix+6]*gtable[:, iy+6]

	ten = bval * (tensors_1 + tensors_2 + tensors_3)
	return ten



## split shape
## NOTE, the bval is a PER GRADIENT scalling
## so here, we have bval=500 and the resulting true bvalue will be 1000 because g1 and g2 have norm 1 and g3 has norm 0
ten_lin = table2ten(grad[:N], bval = 500.)
ten_pla = table2ten(grad[N:], bval = 500.)

## convert mm^2/s to um^2/ms
ten_lin = ten_lin/1000.
ten_pla = ten_pla/1000.



## compute bvalues
b_lin = np.trace(ten_lin, axis1=1, axis2=2)
b_pla = np.trace(ten_pla, axis1=1, axis2=2)


## necessary tensor-product for signal gen
n = np.array([0,0,1]) # tissue orientation
pb_lin = np.outer(n[0]**2, ten_lin[:,0,0]) + np.outer(n[1]**2, ten_lin[:,1,1]) + np.outer(n[2]**2, ten_lin[:,2,2]) + 2*np.outer(n[0]*n[1], ten_lin[:,0,1]) + 2*np.outer(n[0]*n[2], ten_lin[:,0,2]) + 2*np.outer(n[1]*n[2], ten_lin[:,1,2])
pb_pla = np.outer(n[0]**2, ten_pla[:,0,0]) + np.outer(n[1]**2, ten_pla[:,1,1]) + np.outer(n[2]**2, ten_pla[:,2,2]) + 2*np.outer(n[0]*n[1], ten_pla[:,0,1]) + 2*np.outer(n[0]*n[2], ten_pla[:,0,2]) + 2*np.outer(n[1]*n[2], ten_pla[:,1,2])


## non spherical mean signal generation
## V, Vex, Vw := fractions for intra;extra;freewater
def genSig(D1, D2, D3, V, Vex, Vw, pb, b):
    S_i = np.exp(-pb*D1[:,None])
    S_e = np.exp(-pb*D2[:,None] -(b[None,:]-pb)*D3[:,None])
    S_w = np.exp(-b[None,:]*3*np.ones((Vw.shape[0],1)))
    return V[:,None]*S_i + Vex[:,None]*S_e + Vw[:,None]*S_w


## test values
D1 = np.array([0.1])
D2 = np.array([0.6])
D3 = np.array([1.1])
V = np.array([0])
Vex = np.array([1])
Vw = np.array([0])


S_lin = genSig(D1, D2, D3, V, Vex, Vw, pb_lin, b_lin)
S_pla = genSig(D1, D2, D3, V, Vex, Vw, pb_pla, b_pla)


## this should match
signal_lin_numerical = S_lin.mean()
signal_lin_analytical = sm_signal_model_2(b_lin[0], 1., D1, D2, D3, V, Vex, Vw)[0]
signal_pla_numerical = S_pla.mean()
signal_pla_analytical = sm_signal_model_2(b_pla[0], -0.5, D1, D2, D3, V, Vex, Vw)[0]


print('\n')
print('signal linear numerical    {}'.format(signal_lin_numerical))
print('signal linear analytical   {}'.format(signal_lin_analytical))
print('abs diff                   {:.2e}'.format(np.abs(signal_lin_numerical-signal_lin_analytical)))
print('\n')
print('signal linear numerical    {}'.format(signal_pla_numerical))
print('signal linear analytical   {}'.format(signal_pla_analytical))
print('abs diff                   {:.2e}'.format(np.abs(signal_pla_numerical-signal_pla_analytical)))




