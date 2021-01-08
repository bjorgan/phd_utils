#ifndef MNF_C_H_DEFINED
#define MNF_C_H_DEFINED

#include <string>
#include <vector>

enum mnf_err_t {
	MNF_NO_ERR = 0,
	MNF_TRANSFORM_CALCULATION_ERR = -1,
	MNF_STATISTICS_FILE_INCOMPATIBLE = -2,
	MNF_STATISTICS_COVARIANCES_NOT_FOUND = -3,
	MNF_STATISTICS_BANDMEANS_NOT_FOUND = -4,
	MNF_FORWARD_TRANSFORM_CALCULATION_ERR = -5,
	MNF_LU_DECOMPOSITION_ERR = -6,
	MNF_INVERSION_ERR = -7
};

const char *mnf_error_message(mnf_err_t errcode);

/**
 * Image and noise statistics necessary for calculating the MNF transform.
 **/
typedef struct {
	///Image statistics for the noisy image
	void *image_statistics;
	///Statistics for the noise
	void *noise_statistics;
	///For convenience in calculations, number of samples in the noise array
	int num_noise_samples;
	///For convenience in calculations, placeholder array for estimate of noise when doing noise statistics estimation
	float *noise_line;
} mnf_statistics_t;

/**
 * Create mnf_statistics_t structure.
 *
 * \param header Hyperspectral header
 * \return Statistics structure, all fields allocated to the correct sizes
 **/
mnf_statistics_t* mnf_statistics_create(const struct hyspex_header &header);
	
///Postfix for bandmeans file
const std::string MNF_BANDMEANS_FILE_POSTFIX = "_bandmeans.dat";
///Postfix for noise covariance file
const std::string MNF_NOISECOV_FILE_POSTFIX = "_noisecov.dat";
///Postfix for image covariance file
const std::string MNF_IMAGECOV_FILE_POSTFIX = "_imagecov.dat";

/**
 * Check if all MNF statistics files exist (band means, noise and image covariances)
 *
 * \param basefilename Base filename
 * \return True if all files exist, false if all or some are missing
 **/
bool mnf_statistics_exist(std::string basefilename);

/**
 * Read MNF statistics from file.
 *
 * \param filename Input base filename
 * \param header Hyperspectral header
 * \param statistics Return statistics, already allocated using mnf_statistics_create()
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_statistics_from_file(std::string filename, const struct hyspex_header &header, mnf_statistics_t *statistics);

/**
 * Write MNF statistics to file.
 *
 * \param filename Output base filename
 * \param header Hyperspectral header
 * \param statistics MNF statistics
 **/
void mnf_statistics_to_file(std::string filename, const struct hyspex_header &header, const mnf_statistics_t *statistics);

/**
 * Free statistics structure.
 *
 * \param statistics Instance to free
 **/
void mnf_free_statistics(mnf_statistics_t **statistics);

/**
 * Estimate statistics. Will not reinitialize statistics, will (incrementally) update existing statistics with new information. This is done in a numerically stable way.
 *
 * \param statistics Noise and image statistics
 * \param header Hyperspectral header
 * \param img Input hyperspectral image
 **/
void mnf_update_statistics(mnf_statistics_t *statistics, const struct hyspex_header &header, const float *img);

/**
 * MNF transform structure, containing all information necessary for performing forward and inverse MNF transforms.
 **/
typedef struct {
	///Transform from image space to MNF space (transposed, will have to do implicit transpose when applying it to data)
	float *forward_transform;
	///Full inverse transform from MNF space to image space (transposed, will have to do implicit transpose of this matrix when applying it to data)
	float *inverse_transform;
	///Inverse transform with band selection applied (last bands set to zero, so that transform denoises). Not transposed, is already correct.
	float *inverse_transform_with_band_selection;
	///Full denoising transform, yielding denoised spectra from raw, noisy input spectra. Not transposed, already correct.
	float *denoising_transform;
	///Eigenvalues obtained from the eigenvalue problem producing the forward transform
	float *eigenvalues;
	///Band means
	float *means;
	///The matrix R in MNF formulas. Used for selecting the n first bands in the inverse transform, identity matrix with the last m elements set to zero. Used to produce inverse_transform_with_band_selection
	float *band_selection_matrix;
	///Fixed convenience array for operations using the mean
	float *ones_samples;
	///Convenience temporary array for containing the result of a transformation on a hyperspectral line of data
	float *temp_line;
} mnf_transform_t;

/**
 * Create mnf_transform_t structure.
 *
 * \param header Hyperspectral header
 * \param num_bands_in_inverse Number of bands to use in the inverse transform
 * \return Transform structure, with all fields allocated to correct sizes
 **/
mnf_transform_t* mnf_transform_create(const struct hyspex_header &header, int num_bands_in_inverse);

/**
 * Free MNF transform.
 *
 * \param transform MNF instance to free
 **/
void mnf_free_transform(mnf_transform_t **transform);

/**
 * Calculate MNF transform using input statistics.
 *
 * \param header Hyperspectral header
 * \param statistics Image statistics relevant for the MNF transform
 * \param transform Output MNF transform (pre-allocated using mnf_transform_create())
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_transform_calculate(const struct hyspex_header &header, const mnf_statistics_t *statistics, mnf_transform_t *transform);

/**
 * Run MNF forward transform in place on hyperspectral input data, line for line (i.e. each line is treated as a separate matrix which is transformed
 * using the MNF transform).
 *
 * \param transform MNF transform
 * \param header Hyperspectral header
 * \param img Input hyperspectral image. This is overwritten by the result
 **/
void mnf_run_forward(const mnf_transform_t *transform, const struct hyspex_header &header, float *img);

/**
 * Run (band-selected) MNF inverse transform in place on hyperspectral input data.
 *
 * \param transform MNF transform
 * \param header Hyperspectral header
 * \param img MNF transformed image data. This is overwritten by the result
 **/
void mnf_run_inverse(const mnf_transform_t *transform, const struct hyspex_header &header, float *img);

/**
 * Run both forward and inverse transform in place on hyperspectral input data, i.e. will denoise the spectra.
 *
 * \param transform MNF transform
 * \param header Hyperspectral header
 * \param img Hyperspectral image, overwritten by the result
 **/
void mnf_run_full(const mnf_transform_t *transform, const struct hyspex_header &header, float *img);

/**
 * MNF-LBL (incremental line-by-line variant of MNF) algorithm. Will update image/noise statistics
 * with information in the input line of data, recalculate the MNF transform and finally
 * denoise the input data in place.
 *
 * \param transform MNF transform, recalculated using the current statistics
 * \param statistics Image/noise statistics, updated using the new data
 * \param header Hyperspectral header
 * \param line_data Hyperspectral line of data, overwritten by denoised results
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_linebyline(mnf_transform_t *transform, mnf_statistics_t *statistics, const struct hyspex_header &header, float *line_data);

/**
 * Orthonormalize column vectors in the inverse MNF transform, and adjust the forward transform so that the orthonormalized inverse transform can be applied directly.
 * Not used by default in any of the standard routines in this library, but can be used in those cases where it is necessary to do this.
 *
 * \param num_bands_in_inverse Number of bands in inverse (specifies which column vectors to orthornormalize)
 * \param header Hyperspectral image header
 * \param transform Transform to orthornormalize and adjust
 **/
void mnf_orthonormalize_transform(int num_bands_in_inverse, struct hyspex_header header, mnf_transform_t *transform);


#endif
