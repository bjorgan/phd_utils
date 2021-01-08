#ifndef MNF_PRIVATE_H_DEFINED
#define MNF_PRIVATE_H_DEFINED

/**
 * Hyperspectral image covariance and mean.
 **/
typedef struct{
	///Number of pixels summed over
	size_t n;
	///Covariance matrix proto (access using imagestatistics_get_cov(), as some pre-multiplication is needed for true covariance)
	float *C;
	///Band means
	double *means;
	///Fixed convenience array for operations using the mean
	float *ones_samples;
} imagestatistics_t;

imagestatistics_t *image_statistics(mnf_statistics_t *statistics);
imagestatistics_t *noise_statistics(mnf_statistics_t *statistics);

const imagestatistics_t *image_statistics(const mnf_statistics_t *statistics);
const imagestatistics_t *noise_statistics(const mnf_statistics_t *statistics);

/**
 * Calculate forward and inverse MNF transformation matrices.
 *
 * \param bands Number of bands
 * \param imgStats Image statistics
 * \param noiseStats Noise statistics
 * \param forwardTransf Output forward transformation
 * \param inverseTransf Output inverse transformation
 * \param eigvals Output eigenvalues
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_get_transf_matrix(int bands, const imagestatistics_t *imgStats, const imagestatistics_t *noiseStats, float *forwardTransf, float *inverseTransf, float *eigvals);

/**
 * Calculate forward MNF transformation matrix.
 *
 * \param bands Number of bands
 * \param imgCov Image covariance
 * \param noiseCov Noise covariance
 * \param forwardTransf Output forward transformation
 * \param eigvals Eigenvalues
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_calculate_forward_transf_matrix(int bands, const float *imgCov, const float *noiseCov, float *forwardTransf, float *eigvals);

/**
 * Calculate inverse MNF transformation matrix.
 *
 * \param bands Number of bands
 * \param forwardTransf Forward MNF transform
 * \param inverseTransf Outpit inverse transformation
 * \return MNF_NO_ERR on success
 **/
mnf_err_t mnf_calculate_inverse_transf_matrix(int bands, const float *forwardTransf, float *inverseTransf);

/**
 * Estimate noise using the shift difference.
 *
 * \param bands Number of bands
 * \param samples Number of samples
 * \param line Hyperspectral line data
 * \param noise Estimated hyperspectral noise (allocated to at least (num_samples-1)*bands)
 **/
void mnf_estimate_noise_from_line(int bands, int samples, const float *line, float *noise);

/**
 * Remove the mean spectrum from hyperspectral data.
 *
 * \param ones_samples [1, ..., 1], of num_samples length
 * \param means Mean spectrum
 * \param bands Number of bands
 * \param samples Number of samples
 * \param line data, from which the mean spectrum is subtracted
 **/
void mnf_remove_mean_from_line(const float *ones_samples, const float *means, int bands, int samples, float *line);

/**
 * Add the mean spectrum to hyperspectral data.
 *
 * \param ones_samples [1, ..., 1], of num_samples length
 * \param means Mean spectrum
 * \param bands Number of bands
 * \param samples Number of samples
 * \param line data, from which the mean spectrum is subtracted
 **/
void mnf_add_mean_to_line(const float *ones_samples, const float *means, int bands, int samples, float *line);

/**
 * Allocate array fields in imagestatistics_t.
 *
 * \param stats Image statistics to allocate
 * \param samples Number of samples
 * \param bands Number of bands
 **/
void imagestatistics_initialize(imagestatistics_t *stats, int samples, int bands);

/**
 * Deinitialize fields in imagestatistics_t.
 *
 * \param stats Image statistics instance
 **/
void imagestatistics_deinitialize(imagestatistics_t *stats);

/**
 * Update covariance and mean in imagestatistics_t using a line of data in a numerically stable way.
 *
 * \param bands Number of bands
 * \param samples Number of samples
 * \param line_data Hyperspectral line data
 * \param stats Image statistics to update
 **/
void imagestatistics_update_with_line(int bands, int samples, const float *line_data, imagestatistics_t *stats);

/**
 * Obtain mean spectrum from accumulated image statistics.
 *
 * \param stats Image statistics
 * \param bands Number of bands
 * \param means Output mean spectrum
 **/
void imagestatistics_get_means(const imagestatistics_t *stats, int bands, float *means);

/**
 * Obtain covariance matrix from accumulated image statistics.
 *
 * \param stats Image statistics
 * \param bands Number of bands
 * \param cov Output covariance matrix
 **/
void imagestatistics_get_cov(const imagestatistics_t *stats, int bands, float *cov);

/**
 * Read image/noise statistics from file (covariances: only lower row echelon)
 *
 * \param basefilename Base filename
 * \param bands Number of bands
 * \param img_stats Output image statistics
 * \param noise_stats Output noise statistics
 **/
mnf_err_t imagestatistics_read_from_file(std::string basefilename, int bands, imagestatistics_t *img_stats, imagestatistics_t *noise_stats);

/**
 * Write image/noise statistics to file. (Covariances: only lower row echelon)
 *
 * \param basefilename Base filename
 * \param bands Number of bands
 * \param img_stats Image statistics
 * \param noise_stats Noise statistics
 **/
void imagestatistics_write_to_file(std::string basefilename, int bands, const imagestatistics_t *img_stats, const imagestatistics_t *noise_stats);

#endif
