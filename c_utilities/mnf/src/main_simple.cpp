#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "mnf.h"
#include <hyperspectral/readimage.h>
#include "getopt_long_helpers.h"
#include "string_helpers.h"

#define OPT_FORWARD_ONLY 205
#define OPT_INVERSE_ONLY 206

int main(int argc, char *argv[]){
	std::string output_basename;
	int num_bands_in_inverse = 8;
	bool generate_statistics_only = false;
	bool force_statistics_recalculation = false;
	int verbosity = 0;
	bool write_forward_mnf = true;
	bool write_inverse_mnf = true;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"output-basename",		required_argument,	0,	'o'},
		{"statistics-only",		no_argument,		0,	's'},
		{"force-statistics-calculation",no_argument,		0,	'f'},
		{"num-bands-in-inverse",	required_argument,	0,	'b'},
		{"verbose",			no_argument,		0,	'v'},
		{"forward-only", 		no_argument,		0,	OPT_FORWARD_ONLY},
		{"inverse-only",		no_argument,		0,	OPT_INVERSE_ONLY},
		{0, 0, 0, 0}
	};
	char short_options[] = "h:o:s:f:b:v";
	const char *option_descriptions[] = {
		"Show help",
		"Output filename",
		"Write only statistics to file",
		"Force recalculation of image and noise statistics even if the file already exists",
		"Number of bands to use in the inverse transform (default: 8)",
		"Verbose output",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n";

	//parse input arguments
	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;
			case 'o': //output
				output_basename = std::string(optarg);
			break;
			case 's': //generate only statistics
				generate_statistics_only = true;
			break;
			case 'b': //number of bands in inverse
				num_bands_in_inverse = atoi(optarg);
			break;
			case 'f': //force recalculation
				force_statistics_recalculation = true;
			break;
			case 'v':
				verbosity++;
			break;
			case OPT_FORWARD_ONLY:
				write_inverse_mnf = false;
			break;
			case OPT_INVERSE_ONLY:
				write_forward_mnf = false;
			break;
		}
		if (flag == -1){
			break;
		}
	}

	//input filename
	std::string input_filename;
	if (optind >= argc) {
		fprintf(stderr, "Filename missing.\n");
		exit(1);
	} else {
		input_filename = std::string(argv[optind]);
	}

	//output filename
	if (!input_filename.empty() && output_basename.empty()) {
		output_basename = get_basename(input_filename) + "_mnf";
	}
	std::string mnf_forward_filename;
	std::string mnf_inverse_filename;
	if (!output_basename.empty()) {
		mnf_forward_filename = output_basename + "_forwardtransformed";
		mnf_inverse_filename = output_basename + "_inversetransformed";
	}

	//read hyperspectral header
	struct hyspex_header header;
	hyperspectral_err_t errcode;
	errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	//read hyperspectral image
	float *hyperspectral_data = hyperspectral_alloc_float(header);
	errcode = hyperspectral_read_image(input_filename.c_str(), &header, hyperspectral_data);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	//prepare MNF statistics
	mnf_err_t mnf_errcode;
	mnf_statistics_t *statistics = mnf_statistics_create(header);
	mnf_transform_t *transform = mnf_transform_create(header, num_bands_in_inverse);

	if (mnf_statistics_exist(output_basename) && !force_statistics_recalculation) {
		if (verbosity > 0) fprintf(stderr, "Reading MNF statistics from file\n");

		mnf_errcode = mnf_statistics_from_file(output_basename, header, statistics);
		if (mnf_errcode != MNF_NO_ERR) {
			fprintf(stderr, "Error reading statistics files with basename %s: %s\n", output_basename.c_str(), mnf_error_message(mnf_errcode));
			exit(1);
		}
	} else {
		if (verbosity > 0) fprintf(stderr, "Calculating MNF statistics\n");

		mnf_update_statistics(statistics, header, hyperspectral_data);

		if (!output_basename.empty()) {
			if (verbosity > 0) fprintf(stderr, "Writing MNF statistics to file\n");
			mnf_statistics_to_file(output_basename, header, statistics);
		}
	}

	if (!generate_statistics_only) {
		//prepare mnf transform
		mnf_errcode = mnf_transform_calculate(header, statistics, transform);
		if (mnf_errcode != MNF_NO_ERR) {
			fprintf(stderr, "Error calculating MNF transform: %s\n", mnf_error_message(mnf_errcode));
			exit(1);
		}

		hyperspectral_file_t mnf_forward_file, mnf_inverse_file;
		size_t data_size = header.bands*header.samples*header.lines*sizeof(float);

		//run forward MNF
		if (write_forward_mnf) {
			if (verbosity > 0) fprintf(stderr, "Applying forward MNF\n");

			mnf_run_forward(transform, header, hyperspectral_data);
			mnf_forward_file = hyperspectral_open_write_file(mnf_forward_filename, header.bands, header.samples, header.wlens);
			hyperspectral_write_to_file(&mnf_forward_file, data_size, (char*)hyperspectral_data);
			hyperspectral_close_write_file(&mnf_forward_file);
		}

		//run inverse MNF
		if (write_inverse_mnf) {
			if (verbosity > 0) fprintf(stderr, "Applying inverse MNF\n");

			if (write_forward_mnf) {
				mnf_run_inverse(transform, header, hyperspectral_data);
			} else {
				mnf_run_full(transform, header, hyperspectral_data);
			}

			mnf_inverse_file = hyperspectral_open_write_file(mnf_inverse_filename, header.bands, header.samples, header.wlens);
			hyperspectral_write_to_file(&mnf_inverse_file, data_size, (char*)hyperspectral_data);
			hyperspectral_close_write_file(&mnf_inverse_file);
		}

	}

	mnf_free_transform(&transform);
	mnf_free_statistics(&statistics);
}
