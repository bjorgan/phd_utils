#include <unistd.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include "mnf.h"
#include <sys/time.h>
#include <stdio.h>
#include <pthread.h>
#include <sstream>
#include <hyperspectral/readimage.h>
extern "C"
{
#include <lapacke.h>
#include <cblas.h>
}
#include "mnf_private.h"

imagestatistics_t *image_statistics(mnf_statistics_t *statistics)
{
	return (imagestatistics_t*)statistics->image_statistics;
}

imagestatistics_t *noise_statistics(mnf_statistics_t *statistics)
{
	return (imagestatistics_t*)statistics->noise_statistics;
}

const imagestatistics_t *image_statistics(const mnf_statistics_t *statistics)
{
	return (const imagestatistics_t*)statistics->image_statistics;
}

const imagestatistics_t *noise_statistics(const mnf_statistics_t *statistics)
{
	return (const imagestatistics_t*)statistics->noise_statistics;
}

void mnf_update_statistics(mnf_statistics_t *statistics, const struct hyspex_header &header, const float *img)
{
	for (int i=0; i < header.lines; i++) {
		//image statistics
		const float *line_data = img + i*header.samples*header.bands;
		imagestatistics_update_with_line(header.bands, header.samples, line_data, image_statistics(statistics));

		//noise statistics
		mnf_estimate_noise_from_line(header.bands, header.samples, line_data, statistics->noise_line);
		imagestatistics_update_with_line(header.bands, statistics->num_noise_samples, statistics->noise_line, noise_statistics(statistics));
	}
}

bool mnf_statistics_exist(std::string basefilename)
{
	if (access((basefilename + MNF_BANDMEANS_FILE_POSTFIX).c_str(), F_OK) != 0) {
		return false;
	}
	if (access((basefilename + MNF_NOISECOV_FILE_POSTFIX).c_str(), F_OK) != 0) {
		return false;
	}
	if (access((basefilename + MNF_IMAGECOV_FILE_POSTFIX).c_str(), F_OK) != 0) {
		return false;
	}
	return true;
}

mnf_transform_t* mnf_transform_create(const struct hyspex_header &header, int num_bands_in_inverse)
{
	mnf_transform_t transform;
	transform.means = new float[header.bands]();
	transform.forward_transform = new float[header.bands*header.bands]();
	transform.inverse_transform = new float[header.bands*header.bands]();
	transform.eigenvalues = new float[header.bands*header.bands]();
	transform.temp_line = new float[header.samples*header.bands]();
	transform.inverse_transform_with_band_selection = new float[header.bands*header.bands]();
	transform.denoising_transform = new float[header.bands*header.bands]();

	transform.ones_samples = new float[header.samples]();
	for (int i=0; i < header.samples; i++) {
		transform.ones_samples[i] = 1.0f;
	}

	transform.band_selection_matrix = new float[header.bands*header.bands]();
	for (int i=0; i < num_bands_in_inverse; i++) {
		transform.band_selection_matrix[i*header.bands + i] = 1.0f;
	}

	mnf_transform_t *ret_transform = new mnf_transform_t;
	*ret_transform = transform;
	return ret_transform;
}

mnf_statistics_t* mnf_statistics_create(const struct hyspex_header &header)
{
	mnf_statistics_t *ret_stats = new mnf_statistics_t;
	ret_stats->image_statistics = (void*)(new imagestatistics_t);
	ret_stats->noise_statistics = (void*)(new imagestatistics_t);
	ret_stats->num_noise_samples = header.samples-1;
	ret_stats->noise_line = new float[ret_stats->num_noise_samples*header.bands]();
	imagestatistics_initialize(image_statistics(ret_stats), header.samples, header.bands);
	imagestatistics_initialize(noise_statistics(ret_stats), ret_stats->num_noise_samples, header.bands);
	return ret_stats;
}

void mnf_free_statistics(mnf_statistics_t **statistics)
{
	delete [] (*statistics)->noise_line;
	imagestatistics_deinitialize(image_statistics(*statistics));
	imagestatistics_deinitialize(noise_statistics(*statistics));
	delete image_statistics(*statistics);
	delete noise_statistics(*statistics);
	delete *statistics;
}

void mnf_free_transform(mnf_transform_t **transform)
{
	delete [] (*transform)->means;
	delete [] (*transform)->forward_transform;
	delete [] (*transform)->eigenvalues;
	delete [] (*transform)->temp_line;
	delete [] (*transform)->inverse_transform_with_band_selection;
	delete [] (*transform)->denoising_transform;
	delete [] (*transform)->ones_samples;
	delete [] (*transform)->band_selection_matrix;
	delete *transform;
}

void mnf_transform_postmultiply(int bands, const float *band_selection_matrix, const float *forward, const float *inverse, float *inverse_band_selected, float *denoising_transform)
{
	//post-multiply R matrix to get k < m first bands in inverse transformation
	cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, bands, bands, 1.0f, inverse, bands, band_selection_matrix, bands, 0.0f, inverse_band_selected, bands);

	//multiply together the transformation matrices in order to obtain a full denoising matrix
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bands, bands, bands, 1.0f, inverse_band_selected, bands, forward, bands, 0.0f, denoising_transform, bands);
}

mnf_err_t mnf_transform_calculate(const struct hyspex_header &header, const mnf_statistics_t *statistics, mnf_transform_t *transform)
{
	//get means for removing
	imagestatistics_get_means(image_statistics(statistics), header.bands, transform->means);

	//get transfer matrices
	mnf_err_t errcode = mnf_get_transf_matrix(header.bands, image_statistics(statistics), noise_statistics(statistics), transform->forward_transform, transform->inverse_transform, transform->eigenvalues);
	if (errcode != MNF_NO_ERR) {
		return errcode;
	}

	mnf_transform_postmultiply(header.bands, transform->band_selection_matrix, transform->forward_transform, transform->inverse_transform, transform->inverse_transform_with_band_selection, transform->denoising_transform);

	return MNF_NO_ERR;
}

void mnf_run_forward(const mnf_transform_t *transform, const struct hyspex_header &header, float *img)
{
	for (int i=0; i < header.lines; i++) {
		float *line_data = img + i*header.samples*header.bands;

		//remove means
		mnf_remove_mean_from_line(transform->ones_samples, transform->means, header.bands, header.samples, line_data);

		//perform transform of single line
		cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, header.bands, header.samples, header.bands, 1.0f, transform->forward_transform, header.bands, line_data, header.samples, 0.0f, transform->temp_line, header.samples);
		memcpy(line_data, transform->temp_line, sizeof(float)*header.bands*header.samples);
	}
}

void mnf_run_inverse(const mnf_transform_t *transform, const struct hyspex_header &header, float *img)
{
	for (int i=0; i < header.lines; i++) {
		float *line_data = img + i*header.samples*header.bands;

		//perform inverse transform of single line
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, header.bands, header.samples, header.bands, 1.0f, transform->inverse_transform_with_band_selection, header.bands, line_data, header.samples, 0.0f, transform->temp_line, header.samples);

		//add means
		mnf_add_mean_to_line(transform->ones_samples, transform->means, header.bands, header.samples, transform->temp_line);
		memcpy(line_data, transform->temp_line, sizeof(float)*header.samples*header.bands);

	}
}

void mnf_run_full(const mnf_transform_t *transform, const struct hyspex_header &header, float *img)
{
	for (int i=0; i < header.lines; i++) {
		float *line_data = img + i*header.samples*header.bands;

		//remove means
		mnf_remove_mean_from_line(transform->ones_samples, transform->means, header.bands, header.samples, line_data);

		//perform transforms in both direction
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, header.bands, header.samples, header.bands, 1.0f, transform->denoising_transform, header.bands, line_data, header.samples, 0.0f, transform->temp_line, header.samples);

		//add means
		mnf_add_mean_to_line(transform->ones_samples, transform->means, header.bands, header.samples, transform->temp_line);
		memcpy(line_data, transform->temp_line, sizeof(float)*header.samples*header.bands);

	}
}

mnf_err_t mnf_linebyline(mnf_transform_t *transform, mnf_statistics_t *statistics, const struct hyspex_header &full_header, float *line)
{
	struct hyspex_header line_header = full_header;
	line_header.lines = 1; //ensure that all general functions are run on one line only

	//update statistics
	mnf_update_statistics(statistics, line_header, line);

	//calculate transformations
	mnf_err_t errcode = mnf_transform_calculate(line_header, statistics, transform);
	if (errcode != MNF_NO_ERR) {
		return errcode;
	}

	//noise removal
	mnf_run_full(transform, line_header, line);

	return MNF_NO_ERR;
}

mnf_err_t mnf_statistics_from_file(std::string filename, const struct hyspex_header &header, mnf_statistics_t *statistics)
{
	mnf_err_t errcode = imagestatistics_read_from_file(filename, header.bands, image_statistics(statistics), noise_statistics(statistics));
	if (errcode != MNF_NO_ERR) {
		return errcode;
	}
	return MNF_NO_ERR;
}

void mnf_statistics_to_file(std::string filename, const struct hyspex_header &header, const mnf_statistics_t *statistics)
{
	imagestatistics_write_to_file(filename, header.bands, image_statistics(statistics), noise_statistics(statistics));
}

const char *mnf_error_message(mnf_err_t errcode)
{
	switch (errcode) {
		case MNF_LU_DECOMPOSITION_ERR:
			return "LU decomposition failed in calculation of inverse MNF transformation matrix.";
		case MNF_INVERSION_ERR:
			return "Inversion failed in calculation of inverse transformation matrix.";
		case MNF_FORWARD_TRANSFORM_CALCULATION_ERR:
			return "LAPACKE_dsygv failed in calculation of forward MNF transformation matrix.";
	}
	return "Unknown error.";
}

void mnf_orthonormalize_transform(int num_bands_in_inverse, struct hyspex_header header, mnf_transform_t *transform)
{
	//QR decomposition
	float *reflectors = new float[header.bands]();
	int retval = LAPACKE_sgeqrf(LAPACK_COL_MAJOR, header.bands, num_bands_in_inverse, transform->inverse_transform, header.bands, reflectors);
	if (retval != 0) {
		fprintf(stderr, "WARNING: sgeqrf failed.\n");
	}

	//use R to modify the forward transform so that Q^T can be used as the inverse
	float *R = new float[header.bands*header.bands]();
	for (int i=0; i < header.bands; i++) {
		R[i*header.bands + i] = 1.0f;
		for (int j=i; j < num_bands_in_inverse; j++) {
			R[i*header.bands + j] = transform->inverse_transform[j*header.bands + i];
		}
	}

	float *temp = new float[header.bands*header.bands]();
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, header.bands, header.bands, header.bands, 1.0f, R, header.bands, transform->forward_transform, header.bands, 0.0f, temp, header.bands);

	for (int i=0; i < header.bands; i++) {
		for (int j=0; j < header.bands; j++) {
			transform->forward_transform[i*header.bands + j] = temp[j*header.bands + i];
		}
	}
	delete [] temp;
	delete [] R;

	//unpack Q into inverse transform
	retval = LAPACKE_sorgqr(LAPACK_COL_MAJOR, header.bands, num_bands_in_inverse, num_bands_in_inverse, transform->inverse_transform, header.bands, reflectors);
	if (retval != 0) {
		fprintf(stderr, "WARNING: sorgqr failed.\n");
	}

	//prepare the rest of the transforms
	mnf_transform_postmultiply(header.bands, transform->band_selection_matrix, transform->forward_transform, transform->inverse_transform, transform->inverse_transform_with_band_selection, transform->denoising_transform);

	delete [] reflectors;
}

/////////////////////
// Helper routines //
/////////////////////

mnf_err_t mnf_get_transf_matrix(int bands, const imagestatistics_t *imgStats, const imagestatistics_t *noiseStats, float *forwardTransf, float *inverseTransf, float *eigvals)
{
	//get true covariances from supplied statistics
	float *imgCov = new float[bands*bands];
	imagestatistics_get_cov(imgStats, bands, imgCov);

	float *noiseCov = new float[bands*bands];
	imagestatistics_get_cov(noiseStats, bands, noiseCov);

	//estimate forward transformation matrix
	mnf_err_t errcode = mnf_calculate_forward_transf_matrix(bands, imgCov, noiseCov, forwardTransf, eigvals);
	delete [] imgCov;
	delete [] noiseCov;
	if (errcode != MNF_NO_ERR) {
		return errcode;
	}

	//estimate inverse transformation matrix
	errcode = mnf_calculate_inverse_transf_matrix(bands, forwardTransf, inverseTransf);
	if (errcode != MNF_NO_ERR) {
		return errcode;
	}
}

mnf_err_t mnf_calculate_forward_transf_matrix(int bands, const float *imgCov, const float *noiseCov, float *forwardTransf, float *eigvals)
{
	//copy covariances to temp variables, since dsygv will destroy them
	float *noiseCovTemp = new float[bands*bands];
	float *imgCovTemp = new float[bands*bands];
	memcpy(noiseCovTemp, noiseCov, sizeof(float)*bands*bands);
	memcpy(imgCovTemp, imgCov, sizeof(float)*bands*bands);

	//solve eigenvalue problem
	int itype = 1; //solves Ax = lam * Bx
	char jobz = 'V'; //compute both eigenvalues and eigenvectors
	char uplo = 'L'; //using lower triangles of matrices
	int err = LAPACKE_ssygv(LAPACK_ROW_MAJOR, itype, jobz, uplo, bands, noiseCovTemp, bands, imgCovTemp, bands, eigvals);

	//move transfer matrix (stored in noiseCovTemp) to output variable
	memcpy(forwardTransf, noiseCovTemp, sizeof(float)*bands*bands);

	delete [] noiseCovTemp;
	delete [] imgCovTemp;

	if (err != 0) {
		return MNF_FORWARD_TRANSFORM_CALCULATION_ERR;
	}
	return MNF_NO_ERR;

}

mnf_err_t mnf_calculate_inverse_transf_matrix(int bands, const float *forwardTransf, float *inverseTransf)
{
	size_t size = bands*bands*sizeof(float);

	//keep forward transf matrix
	float *forwardTransfTemp = (float*)malloc(size);
	memcpy(forwardTransfTemp, forwardTransf, size);

	//find inverse of forward eigenvector transfer matrix
	int *ipiv = new int[bands];
	int err_sgetrf = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, bands, bands, forwardTransfTemp, bands, ipiv);
	int err_sgetri = LAPACKE_sgetri(LAPACK_ROW_MAJOR, bands, forwardTransfTemp, bands, ipiv);
	delete [] ipiv;

	//copy back
	memcpy(inverseTransf, forwardTransfTemp, size);

	free(forwardTransfTemp);

	if (err_sgetrf != 0) {
		return MNF_LU_DECOMPOSITION_ERR;
	}
	if (err_sgetri) {
		return MNF_INVERSION_ERR;
	}
	return MNF_NO_ERR;
}

void mnf_estimate_noise_from_line(int bands, int samples, const float *line, float *noise)
{
	int num_noise_samples = samples-1;
	for (int i=0; i < bands; i++) {
		for (int j=0; j < num_noise_samples; j++) {
			noise[i*num_noise_samples + j] = line[i*samples + j] - line[i*samples + (j + 1)];
		}
	}
}

void mnf_remove_mean_from_line(const float *ones_samples, const float *means, int bands, int samples, float *line)
{
	cblas_sger(CblasRowMajor, bands, samples, -1.0f, means, 1, ones_samples, 1, line, samples);
}

void mnf_add_mean_to_line(const float *ones_samples, const float *means, int bands, int samples, float *line)
{
	cblas_sger(CblasRowMajor, bands, samples, 1.0f, means, 1, ones_samples, 1, line, samples);
}
