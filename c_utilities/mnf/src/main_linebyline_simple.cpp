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

int main(int argc, char *argv[]){
	std::string output_basename;
	int num_bands_in_inverse = 8;
	int verbosity = 0;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"output-basename",		required_argument,	0,	'o'},
		{"num-bands-in-inverse",	required_argument,	0,	'b'},
		{"verbose",			no_argument,		0,	'v'},
		{0, 0, 0, 0}
	};
	char short_options[] = "h:o:b:v";
	const char *option_descriptions[] = {
		"Show help",
		"Output filename",
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
			case 'b': //number of bands in inverse
				num_bands_in_inverse = atoi(optarg);
			break;
			case 'v':
				verbosity++;
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

	hyperspectral_file_t mnf_file = hyperspectral_open_write_file(output_basename + "_linebyline_denoised", header.bands, header.samples, header.wlens);

	//run mnf-lbl
	for (int i=0; i < header.lines; i++) {
		float *line_data = hyperspectral_data + i*header.samples*header.bands;
		mnf_errcode = mnf_linebyline(transform, statistics, header, line_data);
		if (mnf_errcode != MNF_NO_ERR) {
			fprintf(stderr, "MNF-LBL failed on line %d: %s\n", i, mnf_error_message(mnf_errcode));
			exit(1);
		}

		hyperspectral_write_to_file(&mnf_file, sizeof(float)*header.samples*header.bands, (char*)line_data);
	}
	hyperspectral_close_write_file(&mnf_file);
	mnf_statistics_to_file(output_basename + "_linebyline_statistics", header, statistics);

	mnf_free_transform(&transform);
	mnf_free_statistics(&statistics);
}
