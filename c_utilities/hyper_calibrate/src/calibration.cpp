#include <spectral/spectral.h>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "progress_bar.h"
#include <getopt.h>
#include "calibration.h"
using namespace std;

/**
 * Calculate mean and standard deviation of input data.
 * \param samples Number of samples in data
 * \param data Data array
 * \param output_mean Output mean
 * \param output_std Output standard deviation
 **/
template<typename T>
void get_statistics(int samples, const T *data, float *output_mean, float *output_std){
	double mean = 0;
	double std = 0;
	long n = 0;
	for (int i=0; i < samples; i++){
		n++;
		double delta = data[i] - mean;
		mean = mean + delta/n;
		std = std + delta*(data[i] - mean);
	}

	std = sqrt(std/(n - 1));

	*output_mean = mean;
	*output_std = std;
}

template<typename T>
struct image_subset calibration_find_FOV_refl_standard(enum search_direction direction, struct hyspex_header header, const T *image_data, int verbosity){
	int num_lines = header.lines;
	int num_bands = header.bands;
	int num_samples = header.samples;

	const float MEAN_THRESHOLD = 0.01; //threshold for comparison between subsequent means, defining how much deviation we can have before we can assume the reflectance standard to be finished
	const int MIN_REQUIRED_LINES = 100; //minimum required number of lines in the reflectance standard
	const int NUM_LINES_TO_CHECK_FOR_START = 200;

	int refl_standard_start_line = num_lines;
	int refl_standard_end_line = num_lines;

	//calculate statistics for each line
	int standard_band = num_bands/2;
	float *means = new float[num_lines]();
	float *stds = new float[num_lines]();
	for (int line=0; line < num_lines; line++){
		const T *line_data = image_data + line*num_samples*num_bands;
		get_statistics(num_samples, line_data + standard_band*num_samples, &(means[line]), &(stds[line]));
	}

	//find start line is assumed to always be either start or end of image, 
	//since accounting for the possibility that there might be non-reflectance standard
	//parts of the image at the start or end will make everything too difficult, and
	//will need more advanced search algorithms.
	bool found_start = 0;
	int search_start_line = 1;
	int search_end_line = num_lines - 1;
	int line_increment = 1;
	if (direction == START_FROM_END){
		search_start_line = num_lines-2;
		search_end_line = 0;
		line_increment = -1;
	}
	refl_standard_start_line = search_start_line;
	if (verbosity > 0) fprintf(stderr, "Assuming start line to be %d.\n", refl_standard_start_line);

	for (int line = search_start_line; line_increment*line <= line_increment*search_end_line; line += line_increment){
		int prev_line = line - line_increment;
		double relative_mean = abs((means[line] - means[prev_line]))/(1.0f*means[prev_line]);
		if (relative_mean > MEAN_THRESHOLD){
			refl_standard_end_line = line;
			break;
		}
	}
	if (verbosity > 0) fprintf(stderr, "End line found to be %d\n", refl_standard_end_line);
	delete [] means;
	delete [] stds;

	if (refl_standard_start_line > refl_standard_end_line){
		int temp = refl_standard_end_line;
		refl_standard_end_line = refl_standard_start_line;
		refl_standard_start_line = temp;
	}

	struct image_subset ret_calibration_subset = {0};
	ret_calibration_subset.start_sample = 0;
	ret_calibration_subset.end_sample = num_samples-1;
	ret_calibration_subset.start_band = 0;
	ret_calibration_subset.end_band = num_bands-1;
	ret_calibration_subset.start_line = refl_standard_start_line;
	ret_calibration_subset.end_line = refl_standard_end_line;
	return ret_calibration_subset;
}

struct image_subset calibration_find_FOV_refl_standard(enum search_direction direction, struct hyspex_header header, const char *image_data, int verbosity){
	switch (header.datatype) {
		case HYPERSPECTRAL_DATATYPE_FLOAT:
			return calibration_find_FOV_refl_standard<float>(direction, header, (float*)image_data, verbosity);
		break;

		case HYPERSPECTRAL_DATATYPE_UINT16_T:
			return calibration_find_FOV_refl_standard<uint16_t>(direction, header, (uint16_t*)image_data, verbosity);
		break;
	}
}

hyperspectral_calibration_info_t calibration_create_info(int num_samples, int num_bands, double *calibration_array, struct image_subset refl_standard_positions){
	hyperspectral_calibration_info_t ret_info = {0};
	ret_info.num_samples = num_samples;
	ret_info.num_bands = num_bands;
	ret_info.adjusted = false;
	ret_info.calibration_array = new double[num_samples*num_bands]();
	ret_info.refl_standard_positions = refl_standard_positions;
	memcpy(ret_info.calibration_array, calibration_array, sizeof(double)*num_samples*num_bands);
	return ret_info;
}

void adjust_subset(struct hyspex_header header, struct image_subset &subset)
{
	int num_lines = header.lines;
	int num_bands = header.bands;
	int num_samples = header.samples;

	//adjust subset against image borders
	if (subset.start_sample < 0){
		subset.start_sample = 0;
	}
	if (subset.end_sample >= num_samples){
		subset.end_sample = num_samples-1;
	}
	if (subset.start_line < 0){
		subset.start_line = 0;
	}
	if (subset.end_line >= num_lines){
		subset.end_line = num_lines-1;
	}
}

template<typename T>
hyperspectral_calibration_info_t calibration_integrate(struct hyspex_header header, const T *image_data, struct image_subset subset, int verbosity){
	double *ret_calibration_array = new double[header.samples*header.bands]();

	int num_lines = header.lines;
	int num_bands = header.bands;
	int num_samples = header.samples;

	adjust_subset(header, subset);

	//integrate the image data between start and end line, and let that be the returned calibration spectrum
	long *n = new long[num_bands*num_samples]();
	for (int line = subset.start_line; line <= subset.end_line; line++){
		if (verbosity > 0) print_progress_bar("Integrate reflectance standard.", line - subset.start_line - 1, line - subset.start_line, subset.end_line - subset.start_line - 1);

		for (int band = 0; band < num_bands; band++){
			for (int sample = subset.start_sample; sample <= subset.end_sample; sample++){
				int index = band*num_samples + sample;
				T value = image_data[line*num_bands*num_samples + index];
				n[index]++;
				double delta = value - ret_calibration_array[index];
				ret_calibration_array[index] = ret_calibration_array[index] + delta/n[index];
			}
		}
	}

	for (int band = 0; band < num_bands; band++){
		for (int sample = 0; sample < subset.start_sample; sample++){
			ret_calibration_array[band*num_samples + sample] = ret_calibration_array[band*num_samples + subset.start_sample];
		}
		for (int sample = subset.end_sample + 1; sample < num_samples; sample++){
			ret_calibration_array[band*num_samples + sample] = ret_calibration_array[band*num_samples + subset.end_sample];
		}
	}
	delete [] n;

	bool adjusted = false;
	hyperspectral_calibration_info_t ret_info = calibration_create_info(num_samples, num_bands, ret_calibration_array, subset);

	delete [] ret_calibration_array;

	return ret_info;
}

hyperspectral_calibration_info_t calibration_integrate(struct hyspex_header header, const char *image_data, struct image_subset subset, int verbosity){
	switch (header.datatype) {
		case HYPERSPECTRAL_DATATYPE_FLOAT:
			return calibration_integrate<float>(header, (float*)image_data, subset, verbosity);
		break;

		case HYPERSPECTRAL_DATATYPE_UINT16_T:
			return calibration_integrate<uint16_t>(header, (uint16_t*)image_data, subset, verbosity);
		break;
	}
}

template<typename T>
hyperspectral_calibration_info_t calibration_integrate_all(struct hyspex_header header, const T *image_data, struct image_subset subset, int verbosity){
	int num_lines = header.lines;
	int num_bands = header.bands;
	int num_samples = header.samples;

	adjust_subset(header, subset);

	//integrate the image data
	long *n = new long[num_bands]();
	double *calibration_array = new double[num_bands]();
	for (int line = subset.start_line; line <= subset.end_line; line++){
		if (verbosity > 0) print_progress_bar("Integrate reflectance standard.", line - subset.start_line - 1, line - subset.start_line, subset.end_line - subset.start_line - 1);

		for (int band = 0; band < num_bands; band++){
			for (int sample = subset.start_sample; sample <= subset.end_sample; sample++){
				T value = image_data[line*num_bands*num_samples + band*num_samples + sample];
				n[band]++;
				double delta = value - calibration_array[band];
				calibration_array[band] = calibration_array[band] + delta/n[band];
			}
		}
	}

	delete [] n;

	bool adjusted = false;
	int num_samples_in_calibration_array = 1;
	hyperspectral_calibration_info_t ret_info = calibration_create_info(num_samples_in_calibration_array, num_bands, calibration_array, subset);

	delete [] calibration_array;

	return ret_info;
}

hyperspectral_calibration_info_t calibration_integrate_all(struct hyspex_header header, const char *image_data, struct image_subset subset, int verbosity){
	switch (header.datatype) {
		case HYPERSPECTRAL_DATATYPE_FLOAT:
			return calibration_integrate_all<float>(header, (float*)image_data, subset, verbosity);
		break;

		case HYPERSPECTRAL_DATATYPE_UINT16_T:
			return calibration_integrate_all<uint16_t>(header, (uint16_t*)image_data, subset, verbosity);
		break;
	}
}

void calibration_adjust_band_values(int band, int num_samples, double *calibration_array, double band_specification)
{
	const float CAL_VAL_THRESH = 0.01;
	for (int pix = 0; pix < num_samples; pix++){
		int index = band*num_samples + pix;
		calibration_array[index] /= band_specification*1.0;
		if (calibration_array[index] < CAL_VAL_THRESH){
			calibration_array[index] = CAL_VAL_THRESH;
		}
	}
}

void calibration_adjust(std::string specs_filename, struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info){
	if (!calibration_info->adjusted) {
		double *calibration_array = calibration_info->calibration_array;
		vector<float> wlens = header.wlens;
		int num_bands = calibration_info->num_bands;
		int num_samples = calibration_info->num_samples;

		//read reflectance standard specification file
		spectrum_t specs_spectrum;
		spectral_err_t ret = spectral_read_file(specs_filename.c_str(), &specs_spectrum);
		if (ret != SPECTRAL_NO_ERR){
			printf("Error in reading specification file: %s\n", specs_filename.c_str());
			exit(1);
		}

		//adjust calibration spectrum
		for (int band = 0; band < num_bands; band++){
			float band_specification = 0;
			spectral_get_value(&specs_spectrum, wlens[band], &band_specification);
			calibration_adjust_band_values(band, num_samples, calibration_array, band_specification);
		}

		spectral_free(&specs_spectrum);

		calibration_info->adjusted = true;
	}
}

void calibration_adjust_with_constant_value(double specification_value, struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info){
	if (!calibration_info->adjusted) {
		double *calibration_array = calibration_info->calibration_array;
		vector<float> wlens = header.wlens;
		int num_bands = calibration_info->num_bands;
		int num_samples = calibration_info->num_samples;

		//adjust calibration spectrum
		for (int band = 0; band < num_bands; band++){
			calibration_adjust_band_values(band, num_samples, calibration_array, specification_value);
		}

		calibration_info->adjusted = true;
	}
}

void hyperspectral_calibration_info_free(hyperspectral_calibration_info_t *calibration_info){
	delete [] calibration_info->calibration_array;
}

template<typename T>
void calibration_run_single_line(const hyperspectral_calibration_info_t *calibration_info, struct hyspex_header header, const T *input_data, float *output_data){
	const double *calibration_array = calibration_info->calibration_array;

	for (int j=0; j < header.samples*header.bands; j++){
		int cal_index = j;
		if (calibration_info->num_samples == 1) {
			//calibration array is single-valued
			cal_index = j / header.samples;
		}

		output_data[j] = input_data[j]/calibration_array[cal_index];
	}
}

void calibration_run_single_line(const hyperspectral_calibration_info_t *calibration_info, struct hyspex_header header, const char *input_data, float *output_data){
	switch (header.datatype) {
		case HYPERSPECTRAL_DATATYPE_FLOAT:
			return calibration_run_single_line<float>(calibration_info, header, (float*)input_data, output_data);
		break;

		case HYPERSPECTRAL_DATATYPE_UINT16_T:
			return calibration_run_single_line<uint16_t>(calibration_info, header, (uint16_t*)input_data, output_data);
		break;
	}
}

#ifdef WITH_OPENCV
cv::Mat calibration_generate_rgb_image(struct hyspex_header header, float *rgb_data, struct image_subset subset){
	int num_lines = header.lines;
	int num_bands = header.default_bands.size();
	int num_samples = header.samples;

	std::vector<cv::Mat> channels;
	for (int i=0; (i < num_bands) && (i < 3); i++) {
		channels.push_back(cv::Mat(num_lines, num_samples, CV_32FC1, rgb_data + i*num_samples*num_lines));
	}
	std::reverse(channels.begin(), channels.end());

	cv::Mat rgb_image;
	cv::merge(channels, rgb_image);
	rgb_image = rgb_image.clone()*255;

	//draw calibration subset as rectangle
	cv::Rect cal_subset(subset.start_sample, subset.start_line, subset.end_sample - subset.start_sample, subset.end_line - subset.start_line);
	cv::rectangle(rgb_image, cal_subset, cv::Scalar(0, 0, 255), 10);
	return rgb_image;
}
#endif


template<typename T>
void calibration_crop_refl_standard_in_sample_direction(struct image_subset &image_subset, struct hyspex_header header, const T *image_data, int verbosity){
	int num_lines = image_subset.end_line - image_subset.start_line;
	int num_bands = header.bands;
	int num_samples = header.samples;

	int start_line = image_subset.start_line;
	int end_line = image_subset.end_line;

	int standard_band = num_bands/2;

	float start_line_mean = 0;
	float start_line_stdev;

	float end_line_mean = 0;
	float end_line_stdev = 0;

	const T *start_line_data = image_data + start_line*num_bands*num_samples + standard_band*num_samples;
	const T *end_line_data = image_data + end_line*num_bands*num_samples + standard_band*num_samples;
	get_statistics<T>(num_samples, start_line_data, &start_line_mean, &start_line_stdev);
	get_statistics<T>(num_samples, end_line_data, &end_line_mean, &end_line_stdev);

	bool *start_line_belongs = new bool[num_samples]();
	bool *end_line_belongs = new bool[num_samples]();

	int start_start = -1;
	int start_end = -1;

	int end_start = -1;
	int end_end = -1;

	for (int i=0; i < num_samples; i++) {
		bool start_line_belongs = start_line_data[i] > start_line_mean + start_line_stdev;
		bool end_line_belongs = end_line_data[i] > end_line_mean + end_line_stdev;

		if ((start_start < 0) && (start_line_data[i] > start_line_mean + start_line_stdev)) {
			start_start = i;
		}

		if ((start_end < 0) && (start_line_data[num_samples - i] > start_line_mean + start_line_stdev)) {
			start_end = num_samples - i;
		}

		if ((end_start < 0) && (end_line_data[i] > end_line_mean + end_line_stdev)) {
			end_start = i;
		}

		if ((end_end < 0) && (end_line_data[num_samples - i] > end_line_mean + end_line_stdev)) {
			end_end = num_samples - i;
		}
	}
	image_subset.start_sample = std::max(start_start, end_start);
	image_subset.end_sample = std::min(start_end, end_end);

	if (verbosity > 0) fprintf(stderr, "Start and end samples found to be %d and %d\n", image_subset.start_sample, image_subset.end_sample);


}

void calibration_crop_refl_standard_in_sample_direction(struct image_subset &image_subset, struct hyspex_header header, const char *image_data, int verbosity){
	switch (header.datatype) {
		case HYPERSPECTRAL_DATATYPE_FLOAT:
			calibration_crop_refl_standard_in_sample_direction<float>(image_subset, header, (float*)image_data, verbosity);
		break;

		case HYPERSPECTRAL_DATATYPE_UINT16_T:
			calibration_crop_refl_standard_in_sample_direction<uint16_t>(image_subset, header, (uint16_t*)image_data, verbosity);
		break;
	}
}
