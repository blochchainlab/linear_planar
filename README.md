# Linear-Planar double PGSE sequence processing

This repository contains the preprocessing and fitting scripts for double PGSE data ("Bibek sequence")


### The current preprocessing pipeline uses the following steps

1. Estimate slicewise gaussian sigma (Current version of the pipeline uses an unreleased early version of [Samuel St Jean's Automated characterization of noise distributions in diffusion MRI data](https://www.biorxiv.org/content/10.1101/686436v1))
2. Correct for non-gaussian Bias (second moment method)  
3. Denoising (Mrtrix's implementation of MPPCA)  
4. Gibbs ringing correction (Mrtrix's implementation of the subvoxel shifts method)  
5. Split dataset by b-value and b-tensor shape  
6. Compute susceptibility induced distorsions from AP PA B0 (FSL's topup)  
7. Linear registration (motion correction) and averaging of lowest b-value (FSL's mcflirt)  
8. Get brain mask from moco low b image for eddy current correction (FSL's bet)  
9. Estimate gradient non-linearity distorsion (gradunwarp, requires "secret" coeff.grad file)  
10. Estimate diffusion gradient eddy current induced distosions (FSL's eddy)
11. single interpolation application of (6) (9) (10)  
12. Apply spherical averaging to shell  
13. Compute b0 image from fitting tensor on 2 lowest shells  



### TODOs and Issues
* Because the vector norms in the gradient scheme file are relative, we need to detect the bmax probably from the DICOM header to be able to remove all the hardcoded b-values  
* The one step interpolation of eddy and gradient non-linearity relies on hardcoded software paths  
* the coeff.grad file relies on an hardcoded path  
* Upgrade the gaussian sigma estimation to the [released version](https://github.com/samuelstjean/autodmri)
* Fix and move back to [Cornelius Eichner's version of one step interpolation of eddy and gradient non-linearity](https://github.com/cornelius-eichner/onestep_eddy_nlgc)  





