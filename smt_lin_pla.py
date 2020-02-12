import numpy as np 
from scipy.special import erf
import scipy as sp


# A bit of nomenclature from:
# Lampinen etal., (2019), Searching for the neurite density with diffusion MRI: Challenges for biophysical modeling

# Generic multi-compartement signal formulation
# S = sum(S_k) where S is the signal and S_k are the total for the k compartements

# general form of S_k
# S_k = f_k * S_PD;k * A_T1;k * A_T2;k * A_D;k
# where 
# f_k is the true volume fraction of compartement k
# S_PD;k is the Proton Density of compartement k
# A_T1;k is the T1 attenuation of compartement k
# A_T2;k is the T2 attenuation of compartement k
# A_D;k is the diffusion attenuation of compartement k

# in practice, because we normalize the measured signal by B0
# you either assume equal T1, T2, PD for each compartement and fit "volume fractions"
# or you "do" something like f'_k = f_k * S_PD;k * A_T1;k * A_T2;k
# where f'_k is a signal fraction instead of volume fraction, (i.e. volume fractions weighted by "some exponential" that depend on the underlying tissues parameters)

# Everything here assumes axi-symmetric B-tensor and D-tensor
# so for exemple B = diag(b_para b_perp b_perp) (if this specific B-tensor is "pointing" to [1,0,0])
# so for exemple D = diag(d_para d_perp d_perp) (if this specific D-tensor is "pointing" to [1,0,0])

# because orientation is meaningless for Spherical mean,
# both B-tensor and D-tensor are reparametrize into a "size" parameter and a "shape" parameters (can also be seen as an Isotropic and an Anisotropic parameter)


# D_Iso = Di = (D_para + 2*D_perp) / 3
# D_delta = Dd = (D_para - D_perp) / (D_para + 2*D_perp)

# Same for B-tensor
# b = (normal bval) = trace(B)
# b_delta = bd = {-0.5 if planar; 0 if spherical; 1 if linear}



def main():
	# a few b_delta
	bd_pla = -0.5
	bd_sph =  0.0
	bd_lin =  1.0

	return None



# param conversion utils
def Di_from_param(Dpar, Dperp):
	return (Dpar+2*Dperp)/3.

def Dd_from_param(Dpar, Dperp):
	return (Dpar-Dperp)/float(Dpar+2*Dperp)

def Dpar_from_param(Di, Dd):
	return Di*(1+2*Dd)

def Dperp_from_param(Di, Dd):
	return Di*(1-Dd)




def gfunc(alpha):
	# return sp.sqrt(np.pi/(4*alpha)) * erf(sp.sqrt(alpha))
	tmp = sp.sqrt(np.pi/np.float64(4.*alpha)) * erf(sp.sqrt(alpha))
	return np.where(np.abs(alpha)<1e-16, 1., np.abs(tmp)) # the abs is missing from the paper

# generic signal formula for any B-tensor shape and any D-tensor shape
def sm_signal_generic(b, bd, Di, Dd):
	return np.exp(-b*Di*(1-bd*Dd)) * gfunc(3*b*Di*bd*Dd)

# generic signal formula for any B-tensor shape and any D-tensor shape reparametrized
def sm_signal_generic_reparam(b, bd, Dpar, Dperp):
	Di = Di_from_param(Dpar, Dperp)
	Dd = Dd_from_param(Dpar, Dperp)
	return sm_signal_generic(b, bd, Di, Dd)


# model1 stick(intra) + zepplin(extra) + ball(csf)
def sm_signal_model_1(b, bd, Dpar_in, Dpar_ex, Dperp_ex, f_in, f_ex):
	# stick-like intra axonal compartement with Dperp_ex = 0
	sm_signal_in = sm_signal_generic_reparam(b, bd, Dpar_in, 0)
	# zepplin-like extra axonal compartement
	sm_signal_ex = sm_signal_generic_reparam(b, bd, Dpar_ex, Dperp_ex)
	# ball-like csf axonal compartement with Dpar_ex = Dperp_ex = 3 um^2/ms (free water diffusivity at in-vivo human brain temperature)
	sm_signal_csd = sm_signal_generic_reparam(b, bd, 3., 3.)
	# model1 enforces sum(fraction) = 1 indirectly
	f_csd = 1 - f_in - f_ex
	return f_in*sm_signal_in + f_ex*sm_signal_ex + f_csd*sm_signal_csd

# model2 stick(intra) + zepplin(extra) + ball(csf)
def sm_signal_model_2(b, bd, Dpar_in, Dpar_ex, Dperp_ex, f_in, f_ex, f_csd):
	# stick-like intra axonal compartement with Dperp_ex = 0
	sm_signal_in = sm_signal_generic_reparam(b, bd, Dpar_in, 0)
	# zepplin-like extra axonal compartement
	sm_signal_ex = sm_signal_generic_reparam(b, bd, Dpar_ex, Dperp_ex)
	# ball-like csf axonal compartement with Dpar_ex = Dperp_ex = 3 um^2/ms (free water diffusivity at in-vivo human brain temperature)
	sm_signal_csd = sm_signal_generic_reparam(b, bd, 3., 3.)
	# model2 doesnt enforces sum(fraction) = 1 directly, need to deal with it "later"
	return f_in*sm_signal_in + f_ex*sm_signal_ex + f_csd*sm_signal_csd


