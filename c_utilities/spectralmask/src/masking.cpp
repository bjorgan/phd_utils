#include "masking.h"
#include <spectral/spectral.h>
#include <cmath>
#include <iostream>
#include <string.h>
using namespace std;

masking_err_t masking_init(int num_wlens, std::vector<float*> spectra, masking_t *mask_param, double sam_thresh)
{
	mask_param->num_masking_spectra = spectra.size();
	mask_param->num_bands = num_wlens;
	mask_param->orig_spectra = new float*[spectra.size()];
	mask_param->updated_spectra = new float*[spectra.size()];
	mask_param->sam_thresh = new float[spectra.size()]();
	mask_param->start_band_ind = 0;
	mask_param->end_band_ind = num_wlens - 1;
	mask_param->num_samples_in_spectra = new long[spectra.size()]();

	for (int i=0; i < mask_param->num_masking_spectra; i++){
		mask_param->orig_spectra[i] = new float[num_wlens]();
		mask_param->updated_spectra[i] = new float[num_wlens]();
		memcpy(mask_param->orig_spectra[i], spectra[i], sizeof(float)*num_wlens);
		memcpy(mask_param->updated_spectra[i], spectra[i], sizeof(float)*num_wlens);
		mask_param->sam_thresh[i] = sam_thresh;
	}
	return MASKING_NO_ERR;
}

masking_err_t masking_init(std::vector<float> wlens, spectral_library_t *library, masking_t *mask_param, double sam_thresh)
{
	int num_wlens = wlens.size();

	//generate masking spectra from the spectral library
	std::vector<float*> spectra;
	for (int i=0; i < library->num_spectra; i++) {
		spectra.push_back(new float[num_wlens]());
		for (int j=0; j < num_wlens; j++){
			spectral_get_value(&(library->spectra[i]), wlens[j], &(spectra[i][j]));
		}
	}

	masking_err_t errcode = masking_init(num_wlens, spectra, mask_param, sam_thresh);
	for (int i=0; i < library->num_spectra; i++) {
		delete [] spectra[i];
	}

	return errcode;

}

masking_err_t masking_init(std::vector<float> wlens, std::string spectra_directory, masking_t *mask_param, double sam_thresh)
{
	spectral_library_t library;
	if (spectral_construct_library_from_directory(spectra_directory.c_str(), &library) != SPECTRAL_NO_ERR) {
		return MASKING_LIBRARY_READING_ERR;
	}

	masking_err_t errcode = masking_init(wlens, &library, mask_param, sam_thresh);
	spectral_free_library(&library);
	return errcode;
}

masking_err_t masking_init(std::vector<float> wlens, std::vector<std::string> spectra, masking_t *mask_param, double sam_thresh)
{
	spectral_library_t library;
	char **files = new char*[spectra.size()]();
	for (int i=0; i < spectra.size(); i++) {
		files[i] = new char[spectra[i].length()+1]();
		strcpy(files[i], spectra[i].c_str());
	}

	spectral_err_t spectral_errcode = spectral_construct_library_from_files(spectra.size(), files, &library);
	for (int i=0; i < spectra.size(); i++) {
		delete [] files[i];
	}
	delete [] files;

	if (spectral_errcode != SPECTRAL_NO_ERR) {
		return MASKING_LIBRARY_READING_ERR;
	}

	masking_err_t errcode = masking_init(wlens, &library, mask_param, sam_thresh);
	spectral_free_library(&library);
	return errcode;
}

void masking_free(masking_t *mask_param){
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		delete [] mask_param->orig_spectra[i];
		delete [] mask_param->updated_spectra[i];
	}
	delete [] mask_param->orig_spectra;
	delete [] mask_param->updated_spectra;
	delete [] mask_param->sam_thresh;
	delete [] mask_param->num_samples_in_spectra;
}

void masking_thresh(masking_t *mask_param, int num_samples, float *line_data, mask_thresh_t *ret_thresh){
	//calculate norms of reference spectra
	float *ref_norms_orig = new float[mask_param->num_masking_spectra]();
	float *ref_norms_updated = new float[mask_param->num_masking_spectra]();
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		for (int j=mask_param->start_band_ind; j <= mask_param->end_band_ind; j++){
			ref_norms_orig[i] += mask_param->orig_spectra[i][j]*mask_param->orig_spectra[i][j];
			ref_norms_updated[i] += mask_param->updated_spectra[i][j]*mask_param->updated_spectra[i][j];
		}
		ref_norms_orig[i] = sqrt(ref_norms_orig[i]);
		ref_norms_updated[i] = sqrt(ref_norms_updated[i]);
	}

	for (int j=0; j < num_samples; j++){
		//get pixel band values, calculate norm of pixel spectrum
		float pixel_norm = 0;
		float *pixel_vals = new float[mask_param->num_bands];
		for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
			pixel_vals[i] = line_data[i*num_samples + j];
			pixel_norm += pixel_vals[i]*pixel_vals[i];
		}
		pixel_norm = sqrt(pixel_norm);

		//calculate sam values against all available spectra
		for (int k=0; k < mask_param->num_masking_spectra; k++){
			float samval_orig = 0;
			float samval_updated = 0;
			for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
				samval_orig += pixel_vals[i]*mask_param->orig_spectra[k][i];
				samval_updated += pixel_vals[i]*mask_param->updated_spectra[k][i];
			}
			samval_orig /= pixel_norm*ref_norms_orig[k];
			samval_orig = acos(samval_orig);
			samval_updated /= pixel_norm*ref_norms_updated[k];
			samval_updated = acos(samval_updated);

			//compare against thresholds, save to return array in separate slots
			bool pixel_belong = (samval_orig < mask_param->sam_thresh[k]) || (samval_updated < mask_param->sam_thresh[k]);
			(*ret_thresh)[j][k] = pixel_belong;

			//update the updated spectra with new information if above threshold
			if (pixel_belong){
				long n = mask_param->num_samples_in_spectra[k];
				n++;
				ref_norms_updated[k] = 0;
				for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
					double delta = pixel_vals[i] - mask_param->updated_spectra[k][i];

					//update reference spectrum
					mask_param->updated_spectra[k][i] += delta/(n*1.0);

					//update norm of reference spectrum
					ref_norms_updated[k] += mask_param->updated_spectra[k][i]*mask_param->updated_spectra[k][i];
				}
				ref_norms_updated[k] = sqrt(ref_norms_updated[k]);
				mask_param->num_samples_in_spectra[k] = n;
			}
		}
		delete [] pixel_vals;
	}
	delete [] ref_norms_orig;
	delete [] ref_norms_updated;
}

mask_thresh_t masking_allocate_thresh(const masking_t *mask_param, int num_samples){
	bool **ret_val = new bool*[num_samples];
	for (int i=0; i < num_samples; i++){
		ret_val[i] = new bool[mask_param->num_masking_spectra];
	}
	return ret_val;
}


void masking_free_thresh(mask_thresh_t *mask_thresh, int num_samples){
	for (int i=0; i < num_samples; i++){
		delete [] (*mask_thresh)[i];
	}
	delete [] (*mask_thresh);
}

bool masking_pixel_belongs(const masking_t *mask_param, mask_thresh_t threshed, int sample){
	bool belongs = false;
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		belongs = belongs || threshed[sample][i];
	}
	return belongs;
}

const char *masking_error_message(masking_err_t errcode)
{
	switch (errcode) {
		case MASKING_NO_ERR:
			return "No error.";
		case MASKING_LIBRARY_READING_ERR:
			return "Error in constructing spectral library.";
	}
}
