#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "mnf.h"
#include "../src/mnf_private.h"
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>
#include "getopt_long_helpers.h"
#include "progress_bar.h"
#include "string_helpers.h"

enum opts {
	OPT_NOISE_COVARIANCE = 205,
	OPT_IMAGE_COVARIANCE,
	OPT_BAND_MEAN
};

enum print_opts {
	PRINT_FORWARD,
	PRINT_INVERSE,
	PRINT_REDUCED_INVERSE,
	PRINT_DENOISING,
	PRINT_NOISE_COVARIANCE,
	PRINT_IMAGE_COVARIANCE,
	PRINT_BAND_MEANS,
	PRINT_NONE
};

int main(int argc, char *argv[]){
	int num_bands_in_inverse = 8;
	print_opts print_opt = PRINT_NONE;
	std::string statistics_basename;
	std::string hyperspectral_filename;
	bool apply_orthonormalization = false;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"num-bands-in-inverse",	required_argument,	0,	'b'},
		{"forward-transform",		no_argument,		0,	'f'},
		{"inverse-transform",		no_argument,		0,	'i'},
		{"reduced-inverse-transform",	no_argument,		0,	'r'},
		{"denoising-transform",		no_argument,		0,	'd'},
		{"noise-covariance",		no_argument,		0,	OPT_NOISE_COVARIANCE},
		{"image-covariance",		no_argument,		0,	OPT_IMAGE_COVARIANCE},
		{"band-mean",			no_argument,		0,	OPT_BAND_MEAN},
		{"apply-orthonormalization",	no_argument,		0,	'o'},
		{0, 0, 0, 0}
	};
	char short_options[] = "";
	const char *option_descriptions[] = {
		"Show help",
		"Number of bands to use in the inverse transform (default: 8)",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] STATISTICS_FILE HYPERSPECTRAL_HEADER\n";

	//parse input arguments
	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;

			case 'b': //number of bands in inverse
				num_bands_in_inverse = atoi(optarg);
			break;

			case 'f':
				print_opt = PRINT_FORWARD;
			break;

			case 'i':
				print_opt = PRINT_INVERSE;
			break;

			case 'r':
				print_opt = PRINT_REDUCED_INVERSE;
			break;

			case 'd':
				print_opt = PRINT_DENOISING;
			break;

			case OPT_NOISE_COVARIANCE:
				print_opt = PRINT_NOISE_COVARIANCE;
			break;

			case OPT_IMAGE_COVARIANCE:
				print_opt = PRINT_IMAGE_COVARIANCE;
			break;

			case OPT_BAND_MEAN:
				print_opt = PRINT_BAND_MEANS;
			break;

			case 'o':
				apply_orthonormalization = true;
			break;
		}
		if (flag == -1){
			break;
		}
	}

	//input filename
	if (optind+1 >= argc) {
		fprintf(stderr, "Input statistics filename and header file missing\n");
		exit(1);
	}
	statistics_basename = std::string(argv[optind]);
	hyperspectral_filename = std::string(argv[optind+1]);

	//read hyperspectral header
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header(hyperspectral_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", hyperspectral_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	//prepare MNF statistics
	mnf_err_t mnf_errcode;
	mnf_statistics_t *statistics = mnf_statistics_create(header);
	mnf_errcode = mnf_statistics_from_file(statistics_basename, header, statistics);
	if (mnf_errcode != MNF_NO_ERR) {
		fprintf(stderr, "Error reading statistics files with basename %s: %s\n", statistics_basename.c_str(), mnf_error_message(mnf_errcode));
		exit(1);
	}

	//prepare MNF transform
	mnf_transform_t *transform = mnf_transform_create(header, num_bands_in_inverse);
	mnf_errcode = mnf_transform_calculate(header, statistics, transform);
	if (mnf_errcode != MNF_NO_ERR) {
		fprintf(stderr, "Error calculating MNF transform: %s\n", mnf_error_message(mnf_errcode));
		exit(1);
	}

	//orthornormalize the transform
	if (apply_orthonormalization) {
		mnf_orthonormalize_transform(num_bands_in_inverse, header, transform);
	}

	//obtain statistics arrays for printing
	float *means = new float[header.bands]();
	float *noise_cov = new float[header.bands*header.bands]();
	float *image_cov = new float[header.bands*header.bands]();
	imagestatistics_get_means(image_statistics(statistics), header.bands, means);
	imagestatistics_get_cov(image_statistics(statistics), header.bands, image_cov);
	imagestatistics_get_cov(noise_statistics(statistics), header.bands, noise_cov);

	switch (print_opt) {
		case PRINT_FORWARD:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < header.bands; j++) {
					std::cout << transform->forward_transform[j*header.bands + i] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_INVERSE:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < header.bands; j++) {
					std::cout << transform->inverse_transform[j*header.bands + i] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_REDUCED_INVERSE:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < header.bands; j++) {
					std::cout << transform->inverse_transform_with_band_selection[i*header.bands + j] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_DENOISING:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < header.bands; j++) {
					std::cout << transform->denoising_transform[i*header.bands + j] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_NOISE_COVARIANCE:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < i; j++) {
					std::cout << noise_cov[i*header.bands + j] << " ";
				}

				for (int j=i; j < header.bands; j++) {
					std::cout << noise_cov[j*header.bands + i] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_IMAGE_COVARIANCE:
			for (int i=0; i < header.bands; i++) {
				for (int j=0; j < i; j++) {
					std::cout << image_cov[i*header.bands + j] << " ";
				}

				for (int j=i; j < header.bands; j++) {
					std::cout << image_cov[j*header.bands + i] << " ";
				}
				std::cout << std::endl;
			}
		break;

		case PRINT_BAND_MEANS:
			for (int i=0; i < header.bands; i++) {
				std::cout << header.wlens[i] << " " << means[i] << std::endl;
			}
		break;

		default:
			fprintf(stderr, "Nothing to print.\n");
		break;
	}

	delete [] means;
	delete [] noise_cov;
	delete [] image_cov;
	mnf_free_transform(&transform);
	mnf_free_statistics(&statistics);
}
