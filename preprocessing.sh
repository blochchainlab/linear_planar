# This script assumes:
    # Scheme encoding file: $FOLDER'G1_G2_G3.txt'
    # Modified topup option file: $FOLDER'b02b0_2mm.cnf' # TODO: move this outside of the data folders
    # a folder $FOLDER'NII/' containing 
        # dmri data: dwi.nii
        # eddy config file: eddy_config.txt
        # eddy indicies file: eddy_indices.txt
        # AP b0: AP.nii
        # PA b0: PA.nii 
        # Phase Encoding table AP: pe_table_AP.txt
        # Phase Encoding table PA: pe_table_PA.txt

# relies on a few environement: fsl5, mrtrix3, ants, cuda8

# for eddy_openmp
OMP_NUM_THREADS=6

# CHANGE ME
FOLDER=/path/to/root/
GRADFILE=coeffConnectom.grad

# path to https://github.com/cornelius-eichner/onestep_eddy_nlgc
ONESHOTPATH=/home/raid2/paquette/tools/onestep_eddy_nlgc/

# Make folder where all the (very numerous) sub-step will live
mkdir -p $FOLDER'PREPRO'

# raw sigma estimation (with or without noisemap)
get_distribution $FOLDER'NII/dwi.nii' $FOLDER'PREPRO/dwi_SAM_sigma.nii' $FOLDER'PREPRO/dwi_SAM_N.nii' $FOLDER'PREPRO/dwi_SAM_mask.nii' -a 2 -m 'maxlk' --ncores 6
# get_distribution $FOLDER'NII/dwi.nii' $FOLDER'PREPRO/dwi_SAM_sigma.nii' $FOLDER'PREPRO/dwi_SAM_N.nii' $FOLDER'PREPRO/dwi_SAM_mask.nii' -a 2 -m 'maxlk' --ncores 6 --noise_maps

# rician bias correction
debias.py $FOLDER'NII/dwi.nii' $FOLDER'PREPRO/dwi_SAM_N.nii' $FOLDER'PREPRO/dwi_SAM_sigma.nii' $FOLDER'PREPRO/dwi_deb.nii' $FOLDER'PREPRO/dwi_deb_nan.nii' $FOLDER'PREPRO/dwi_deb_nanMean.nii'

# denoising
dwidenoise $FOLDER'PREPRO/dwi_deb.nii' $FOLDER'PREPRO/dwi_deb_den.nii' -noise $FOLDER'PREPRO/dwi_deb_mppca_sigma.nii'
mrcalc $FOLDER'PREPRO/dwi_deb.nii' $FOLDER'PREPRO/dwi_deb_den.nii' -subtract $FOLDER'PREPRO/dwi_deb_den_Residual.nii'
mrabs $FOLDER'PREPRO/dwi_deb_den_Residual.nii' $FOLDER'PREPRO/dwi_deb_den_ResidualAbs.nii'

# gibbs ringing correction
mrdegibbs $FOLDER'PREPRO/dwi_deb_den.nii' $FOLDER'PREPRO/dwi_deb_den_gib.nii' -axes 0,1
mrcalc $FOLDER'PREPRO/dwi_deb_den.nii' $FOLDER'PREPRO/dwi_deb_den_gib.nii' -subtract $FOLDER'PREPRO/dwi_deb_den_gib_Residual.nii'
mrabs $FOLDER'PREPRO/dwi_deb_den_gib_Residual.nii' $FOLDER'PREPRO/dwi_deb_den_gib_ResidualAbs.nii'

# split data by b-value b-tensor shape for shell by shell eddy
# TODO: make automatic detection of b-values
split_shells_and_types.py $FOLDER'PREPRO/dwi_deb_den_gib.nii' 50,1000,2000,4000 $FOLDER'G1_G2_G3.txt' 4000 $FOLDER'NII/eddy_config.txt' $FOLDER'NII/eddy_indices.txt' $FOLDER'PREPRO/dwi_deb_den_gib' $FOLDER'PREPRO/bvecs' $FOLDER'PREPRO/bvals' $FOLDER'PREPRO/indices' $FOLDER'PREPRO/positions'

# topup 
mrcat -axis 3 $FOLDER'NII/AP.nii' $FOLDER'NII/PA.nii' $FOLDER'NII/AP_PA.nii'
cat $FOLDER'NII/pe_table_AP.txt' $FOLDER'NII/pe_table_PA.txt' > $FOLDER'NII/pe_table_AP_PA.txt'
topup --imain=$FOLDER'NII/AP_PA.nii' --datain=$FOLDER'NII/pe_table_AP_PA.txt' --out=$FOLDER'PREPRO/topup_' --iout=$FOLDER'PREPRO/topup.nii' --config=$FOLDER'b02b0_2mm.cnf'

# Linear registration of spreaded out b50
fsl5.0-mcflirt -in $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin.nii' -out $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_mc.nii.gz' -stages 4 -refvol 0 

# average motion corrected b50
mrmath $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_mc.nii.gz' mean $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_mc_mean.nii.gz' -axis 3

# make quick brain mask from mean b50
# This bet should definitely be manually QCd
bet $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_mc_mean.nii.gz' $FOLDER'PREPRO/bet_init_b50.nii.gz' -m -n -f 0.1
maskfilter $FOLDER'PREPRO/bet_init_b50_mask.nii.gz' dilate $FOLDER'PREPRO/bet_init_b50_dil.nii.gz'

# fsl EDDY shell by shell with GNLC single interpolation
# TODO: not hardcoding b-values
for BVAL in 1000 2000 4000;
do
for SHELLTYPE in lin pla;
do 
echo 'EDDY for b = '$BVAL'  '$SHELLTYPE;
# pre-append mean b50 and make tmp indices, bvec and bval accordingly 
mrcat -axis 3 $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_mc_mean.nii.gz' $FOLDER'PREPRO/dwi_deb_den_gib__b'$BVAL'_'$SHELLTYPE'.nii' $FOLDER'PREPRO/preeddy__b'$BVAL'_'$SHELLTYPE'.nii';
(echo 1; cat $FOLDER'PREPRO/indices__b'$BVAL'_'$SHELLTYPE'.txt') > $FOLDER'PREPRO/preeddy_indices__b'$BVAL'_'$SHELLTYPE'.txt';
# depends if bvec are rows or columns
# sed -e 's/^/0.0 /' $FOLDER'PREPRO/bvecs__b'$BVAL'_'$SHELLTYPE'.txt' > $FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt'
(echo 0.0 0.0 0.0; cat $FOLDER'PREPRO/bvecs__b'$BVAL'_'$SHELLTYPE'.txt') > $FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt';
(echo 0.0; cat $FOLDER'PREPRO/bvals__b'$BVAL'_'$SHELLTYPE'.txt') > $FOLDER'PREPRO/preeddy_bvals__b'$BVAL'_'$SHELLTYPE'.txt';
# uncomment to use eddy openmp without gradient non-linearity correction
# /home/raid2/paquette/tools/eddy_openmp --imain=$FOLDER'PREPRO/preeddy__b'$BVAL'_'$SHELLTYPE'.nii' --mask=$FOLDER'PREPRO/bet_init_b50_dil.nii.gz' --index=$FOLDER'PREPRO/preeddy_indices__b'$BVAL'_'$SHELLTYPE'.txt' --acqp=$FOLDER'NII/pe_table_AP_PA.txt' --bvecs=$FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt' --bvals=$FOLDER'PREPRO/preeddy_bvals__b'$BVAL'_'$SHELLTYPE'.txt' --out=$FOLDER'PREPRO/eddy_b'$BVAL'_'$SHELLTYPE'' --topup=$FOLDER'PREPRO/topup_' --cnr_maps --data_is_shelled; 
# uncomment to use eddy cuda without gradient non-linearity correction
# /home/raid2/paquette/tools/eddy_cuda8.0 --imain=$FOLDER'PREPRO/preeddy__b'$BVAL'_'$SHELLTYPE'.nii' --mask=$FOLDER'PREPRO/bet_init_b50_dil.nii.gz' --index=$FOLDER'PREPRO/preeddy_indices__b'$BVAL'_'$SHELLTYPE'.txt' --acqp=$FOLDER'NII/pe_table_AP_PA.txt' --bvecs=$FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt' --bvals=$FOLDER'PREPRO/preeddy_bvals__b'$BVAL'_'$SHELLTYPE'.txt' --out=$FOLDER'PREPRO/eddy_b'$BVAL'_'$SHELLTYPE'' --topup=$FOLDER'PREPRO/topup_' --cnr_maps --data_is_shelled; 
# make sure it MULTIPLY the jacobian
python $ONESHOTPATH"eddy_nlgc_mp.py" --in $FOLDER'PREPRO/preeddy__b'$BVAL'_'$SHELLTYPE'.nii' --bvec $FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt' --bval $FOLDER'PREPRO/preeddy_bvals__b'$BVAL'_'$SHELLTYPE'.txt' --mask $FOLDER'PREPRO/bet_init_b50_dil.nii.gz' --acqp $FOLDER'NII/pe_table_AP_PA.txt' --index $FOLDER'PREPRO/preeddy_indices__b'$BVAL'_'$SHELLTYPE'.txt' --topup $FOLDER'PREPRO/topup_' --out $FOLDER'PREPRO/onestep__b'$BVAL'_'$SHELLTYPE'.nii'
# cleanup to keep everything per shell
mv $FOLDER'PREPRO/eddy_tmp.eddy_post_eddy_shell_PE_translation_parameters' $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_post_eddy_shell_PE_translation_parameters'  
mv $FOLDER'PREPRO/eddy_tmp.eddy_post_eddy_shell_alignment_parameters'      $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_post_eddy_shell_alignment_parameters'       
mv $FOLDER'PREPRO/eddy_tmp.eddy_outlier_report'                            $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_outlier_report'                             
mv $FOLDER'PREPRO/eddy_tmp.eddy_outlier_map'                               $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_outlier_map'                                                        
mv $FOLDER'PREPRO/eddy_tmp.eddy_outlier_n_stdev_map'                       $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_outlier_n_stdev_map'                                                        
mv $FOLDER'PREPRO/eddy_tmp.eddy_outlier_n_sqr_stdev_map'                   $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_outlier_n_sqr_stdev_map'                                                
mv $FOLDER'PREPRO/eddy_tmp.eddy_outlier_free_data.nii.gz'                  $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_outlier_free_data.nii.gz'                                               
mv $FOLDER'PREPRO/eddy_tmp.eddy_parameters'                                $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_parameters'                                                 
mv $FOLDER'PREPRO/eddy_tmp.nii.gz'                                         $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.nii.gz'                                                  
mv $FOLDER'PREPRO/eddy_tmp.eddy_rotated_bvecs'                             $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_rotated_bvecs'                                                  
mv $FOLDER'PREPRO/eddy_tmp.eddy_movement_rms'                              $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_movement_rms'                                                   
mv $FOLDER'PREPRO/eddy_tmp.eddy_restricted_movement_rms'                   $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_restricted_movement_rms'                                                    
mv $FOLDER'PREPRO/eddy_tmp.eddy_cnr_maps.nii.gz'                           $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_eddy_tmp.eddy_cnr_maps.nii.gz'                                
rm $FOLDER'PREPRO/single.nii.gz' -f
mv $FOLDER'PREPRO/nlgc_warp.nii.gz'                                        $FOLDER'PREPRO/nlgc_warp_b'$BVAL'_'$SHELLTYPE'.nii.gz'
mv $FOLDER'PREPRO/nlgc.nii.gz'                                             $FOLDER'PREPRO/nlgc_b'$BVAL'_'$SHELLTYPE'.nii.gz'
mv $FOLDER'PREPRO/dfield/'                                                 $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_dfield/'
mv $FOLDER'PREPRO/comb_warp/'                                              $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_comb_warp/'
mv $FOLDER'PREPRO/split_data/'                                             $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_split_data/'
mv $FOLDER'PREPRO/warped_data/'                                            $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_warped_data/'
mv $FOLDER'PREPRO/jac_det_data/'                                           $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_jac_det_data/'
mv $FOLDER'PREPRO/corr_data/'                                              $FOLDER'PREPRO/b'$BVAL'_'$SHELLTYPE'_corr_data/'
dwiextract -no_bzero -fslgrad $FOLDER'PREPRO/preeddy_bvecs__b'$BVAL'_'$SHELLTYPE'.txt' $FOLDER'PREPRO/preeddy_bvals__b'$BVAL'_'$SHELLTYPE'.txt' $FOLDER'PREPRO/onestep__b'$BVAL'_'$SHELLTYPE'.nii.gz' $FOLDER'PREPRO/v1_eddy_b'$BVAL'_'$SHELLTYPE'.nii.gz';
done;
done

# compute spherical mean
# TODO: not hardcoding b-values
for BVAL in 1000 2000 4000;
do
for SHELLTYPE in lin pla;
do 
echo 'SM for b = '$BVAL'  '$SHELLTYPE;
mean_std_4th.py $FOLDER'PREPRO/v1_eddy_b'$BVAL'_'$SHELLTYPE'.nii.gz' $FOLDER'PREPRO/v1_eddy_b'$BVAL'_'$SHELLTYPE'_SM.nii.gz' $FOLDER'PREPRO/v1_eddy_b'$BVAL'_'$SHELLTYPE'_SM_std.nii.gz'
done;
done

# making a S0 image from b50 lin and b1000 lin
# 1) apply topup to b50
applytopup -i $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin.nii' -a $FOLDER'NII/pe_table_AP_PA.txt' -x 1 -t $FOLDER'PREPRO/topup_' -o $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top.nii.gz' -m jac
# 2) run grdient non linear correction
gradient_unwarp.py $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top.nii.gz' $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc_nojac.nii' siemens -g $GRADFILE -n --outfilewarp $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_warp.nii.gz'
python $ONESHOTPATH"calc_jacobian.py" --in $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_warp.nii.gz' --out $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_warp_jac.nii.gz'
fslmaths $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc_nojac.nii' -mul $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_warp_jac.nii.gz' $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc.nii'
# 3) run mc-flirt on topuped b50
fsl5.0-mcflirt -in $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc.nii.gz' -out $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc_mc.nii' -stages 4 -refvol 0
# 4) cat with topup+eddy corrected b1000
mrcat -axis 3 $FOLDER'PREPRO/dwi_deb_den_gib__b50_lin_top_nlgc_mc.nii.gz' $FOLDER'PREPRO/v1_eddy_b1000_lin.nii.gz' $FOLDER'PREPRO/dwi_for_S0.nii'
(cat $FOLDER'PREPRO/bvecs__b50_lin.txt'; cat $FOLDER'PREPRO/bvecs__b1000_lin.txt') > $FOLDER'PREPRO/bvecs_for_S0.txt'
(cat $FOLDER'PREPRO/bvals__b50_lin.txt'; cat $FOLDER'PREPRO/bvals__b1000_lin.txt') > $FOLDER'PREPRO/bvals_for_S0.txt'
# 5) fit tensor to get S0
dwi2tensor $FOLDER'PREPRO/dwi_for_S0.nii' $FOLDER'PREPRO/dti_for_S0.nii' -b0 $FOLDER'PREPRO/computed_S0.nii' -fslgrad $FOLDER'PREPRO/bvecs_for_S0.txt' $FOLDER'PREPRO/bvals_for_S0.txt' -mask $FOLDER'PREPRO/bet_init_b50_dil.nii.gz'




