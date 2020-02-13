import numpy as np
import nibabel as nib


def main(fdata, fNs, fsigmas, fcorr_data, fnan_map, fnan_mean):
    data_img = nib.load(fdata)
    data = data_img.get_data()

    Ns = nib.load(fNs).get_data()
    sigmas = nib.load(fsigmas).get_data()

    corr_data = np.sqrt(data**2 - 2*Ns[...,None]*(sigmas[...,None]**2))

    nan_map = np.zeros_like(corr_data, dtype=np.bool)
    nan_idx = np.where(np.isnan(corr_data))
    nan_map[nan_idx] = True

    corr_data[nan_idx] = 0

    nib.Nifti1Image(corr_data, data_img.affine, data_img.header).to_filename(fcorr_data)
    nib.Nifti1Image(nan_map, data_img.affine, data_img.header).to_filename(fnan_map)
    nib.Nifti1Image(nan_map.mean(3), data_img.affine, data_img.header).to_filename(fnan_mean)


if __name__ == "__main__":
    import sys
    print('usage: python debias.py input_data input_Ns input_sigmas output_corr_data output_nan_map nana_mean_map')
    main(*sys.argv[1:])
