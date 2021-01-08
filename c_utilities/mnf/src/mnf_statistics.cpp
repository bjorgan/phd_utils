#include "mnf.h"
#include "mnf_private.h"
#include <fstream>
#include <string.h> 

extern "C" {
#include <cblas.h>
}
#include <sstream>
#include <limits>

void imagestatistics_get_means(const imagestatistics_t *stats, int numBands, float *means)
{
	for (int i=0; i < numBands; i++) {
		means[i] = stats->means[i];
	}
}

void imagestatistics_get_cov(const imagestatistics_t *stats, int numBands, float *cov)
{
	for (int i=0; i < numBands*numBands; i++) {
		cov[i] = stats->C[i]/(stats->n*1.0f);
	}
}

mnf_err_t imagestatistics_read_from_file(std::string basefilename, int bands, imagestatistics_t *imgStats, imagestatistics_t *noiseStats)
{
	std::string imgCovStr = basefilename + MNF_IMAGECOV_FILE_POSTFIX;
	std::string noiseCovStr = basefilename + MNF_NOISECOV_FILE_POSTFIX;

	std::ifstream imgCovFile;
	imgCovFile.open(imgCovStr.c_str());

	std::ifstream noiseCovFile;
	noiseCovFile.open(noiseCovStr.c_str());

	if (imgCovFile.fail() || imgCovFile.fail()) {
		return MNF_STATISTICS_COVARIANCES_NOT_FOUND;
	}

	//copy to container
	for (int i=0; i < bands; i++) {
		std::string imgLine, noiseLine;
		getline(imgCovFile, imgLine);
		getline(noiseCovFile, noiseLine);
		std::stringstream sI(imgLine);
		std::stringstream sN(noiseLine);
		for (int j=0; j < bands; j++) {
			float val;
			sI >> val;
			imgStats->C[i*bands + j] = val;
			sN >> val;
			noiseStats->C[i*bands + j] = val;
		}
	}

	imgStats->n = 1;
	noiseStats->n = 1;

	noiseCovFile.close();
	imgCovFile.close();


	//read band means from file
	std::string meanStr = basefilename + MNF_BANDMEANS_FILE_POSTFIX;
	std::ifstream meanFile;
	meanFile.open(meanStr.c_str());

	if (meanFile.fail()) {
		return MNF_STATISTICS_BANDMEANS_NOT_FOUND;
	}

	//copy to container
	for (int i=0; i < bands; i++) {
		float val = 0;
		meanFile >> val;
		imgStats->means[i] = val;
	}

	meanFile.close();

	return MNF_NO_ERR;
}

void imagestatistics_write_cov_to_file(const imagestatistics_t *imgStats, int bands, std::string filename)
{
	std::ofstream file;
	file.open(filename.c_str());
	float *cov = new float[bands*bands];
	imagestatistics_get_cov(imgStats, bands, cov);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	for (int i=0; i < bands; i++) {
		for (int j=0; j < bands; j++) {
			file << cov[i*bands + j] << " ";
		}
		file << std::endl;
	}
	file.close();
	delete [] cov;
}

void imagestatistics_write_mean_to_file(const imagestatistics_t *imgStats, int bands, std::string filename)
{
	float *means = new float[bands];
	imagestatistics_get_means(imgStats, bands, means);
	std::ofstream file;
	file.open(filename.c_str());
	for (int i=0; i < bands; i++) {
		file << means[i] << " ";
	}
	file << std::endl;
	file.close();
	delete [] means;
}

void imagestatistics_write_to_file(std::string basefilename, int bands, const imagestatistics_t *imgStats, const imagestatistics_t *noiseStats)
{
	imagestatistics_write_cov_to_file(imgStats, bands, basefilename + MNF_IMAGECOV_FILE_POSTFIX);
	imagestatistics_write_cov_to_file(noiseStats, bands, basefilename + MNF_NOISECOV_FILE_POSTFIX);
	imagestatistics_write_mean_to_file(imgStats, bands, basefilename + MNF_BANDMEANS_FILE_POSTFIX);
}

void imagestatistics_initialize(imagestatistics_t *stats, int samples, int bands)
{
	stats->n = 0;
	stats->C = new float[bands*bands];
	for (int i=0; i < bands*bands; i++) {
		stats->C[i] = 0.0f;
	}

	stats->means = new double[bands];
	for (int i=0; i < bands; i++) {
		stats->means[i] = 0.0f;
	}

	stats->ones_samples = new float[samples]();
	for (int i=0; i < samples; i++) {
		stats->ones_samples[i] = 1.0f;
	}
}

void imagestatistics_deinitialize(imagestatistics_t *stats)
{
	delete [] stats->C;
	delete [] stats->means;
}

void imagestatistics_update_with_line(int bands, int samples, const float *bilData, imagestatistics_t *stats)
{
	stats->n += samples;
	float *C = stats->C;
	double *means = stats->means;

	//copy line to temporary variable
	float *tempLine = new float[samples*bands];
	memcpy(tempLine, bilData, sizeof(float)*samples*bands);

	//estimate mean of single line
	float *meanTemp = new float[bands];
	cblas_sgemv(CblasRowMajor, CblasNoTrans, bands, samples, 1.0f/(1.0f*samples), tempLine, samples, stats->ones_samples, 1, 0.0f, meanTemp, 1);

	//subtract mean from line
	mnf_remove_mean_from_line(stats->ones_samples, meanTemp, bands, samples, tempLine);

	//calculate covariance for current line and add to accumulated covariance
	cblas_ssyrk(CblasRowMajor, CblasLower, CblasNoTrans, bands, samples, 1.0f, tempLine, samples, 1.0f, C, bands);

	//find the difference between the mean of the current line and the old mean
	//update total means
	float *meanDiff = new float[bands];
	for (int i=0; i < bands; i++) {
		meanDiff[i] = meanTemp[i] - means[i];
		means[i] = means[i] + 1.0*samples*(meanTemp[i] - means[i])/(1.0*stats->n);
	}

	//update to actual covariance
	cblas_ssyr(CblasRowMajor, CblasLower, bands, samples*(stats->n - samples)/(stats->n), meanDiff, 1, C, bands);

	delete [] meanDiff;
	delete [] meanTemp;
	delete [] tempLine;
}
