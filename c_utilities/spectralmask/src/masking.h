#ifndef MASKING_H_DEFINED
#define MASKING_H_DEFINED

#include <vector>
#include <string>

/**
 * Masking parameters. Reference spectra and so on.
 **/
typedef struct{
	/// Number of reference spectra
	int num_masking_spectra;
	/// Number of bands
	int num_bands;
	/// Original reference spectra, as input into the initializator
	float **orig_spectra;
	/// Spectra that are updated with new information as the image is segmented as skin
	float **updated_spectra;
	/// Number of samples used in updated_spectra
	long *num_samples_in_spectra;
	/// Threshold values for SAM
	float *sam_thresh;
	/// Start band for SAM calculations
	int start_band_ind;
	/// End band for SAM calculations
	int end_band_ind;
} masking_t;

/**
 * Masking error values.
 **/
enum masking_err_t{
	MASKING_NO_ERR = 0,
	MASKING_LIBRARY_READING_ERR = -1
};

/**
 * Masking error messages.
 **/
const char *masking_error_message(masking_err_t errcode);

/**
 * Internal datatype for controlling segmentations. Samples along the first table direction, ref. spectrum number along the next.
 **/
typedef bool** mask_thresh_t;

/**
 * Allocate mask_thresh_t object.
 * \param mask_param Masking parameters
 * \param num_samples Number of samples in image
 **/
mask_thresh_t masking_allocate_thresh(const masking_t *mask_param, int num_samples);

/**
 * Free mask_thresh_t object. Everything will be freed, no need to do anything from outside.
 **/
void masking_free_thresh(mask_thresh_t *mask_thresh_t, int num_samples);

/**
 * Check whether specified pixel belongs according to the thresholded SAM values.
 * \param mask_param Masking parameters
 * \param threshed Thresholded values obtained from masking_thresh()
 * \param sample Sample coordinate along image
 * \return true if pixel belongs to the segmented image
 **/
bool masking_pixel_belongs(const masking_t *mask_param, mask_thresh_t threshed, int sample);

#define SAM_THRESH_DEFAULT 0.3

/**
 * Initialize masking parameters.
 * \param num_wlens Number of bands in image to segment
 * \param spectra List of spectra
 * \param mask_param Output masking parameters
 * \param sam_thresh SAM threshold
 * \return MASKING_NO_ERR on success
 **/
masking_err_t masking_init(int num_wlens, std::vector<float*> spectra, masking_t *mask_param, double sam_thresh);

/**
 * Initialize masking parameters.
 * \param wlens Wavelengths
 * \param spectra_directory Path to directory containing masking spectra
 * \param mask_param Output masking parameters
 * \param sam_thresh SAM threshold
 * \return MASKING_NO_ERR on success
 **/
masking_err_t masking_init(std::vector<float> wlens, std::string spectra_directory, masking_t *mask_param, double sam_thresh);

/**
 * Initialize masking parameters.
 * \param wlens Wavelengths
 * \param spectra List over filenames for masking spectra
 * \param mask_param Output masking parameters
 * \param sam_thresh SAM threshold
 * \return MASKING_NO_ERR on success
 **/
masking_err_t masking_init(std::vector<float> wlens, std::vector<std::string> spectra, masking_t *mask_param, double sam_thresh);

/**
 * Do masking thresholding according to parameter specifications and update reference spectra according to segmented parts.
 * \param mask_param Masking parameters
 * \param num_samples Number of samples in image
 * \param line_data Input hyperspectral data
 * \param ret_thresh Return segmented values.
 **/
void masking_thresh(masking_t *mask_param, int num_samples, float *line_data, mask_thresh_t *ret_thresh);

/**
 * Free memory associated with masking parameters.
 **/
void masking_free(masking_t *mask_param);

#endif
