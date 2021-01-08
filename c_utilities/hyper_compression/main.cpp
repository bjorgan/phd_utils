#include "compression.h"
#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>
#include "getopt_long_helpers.h"
#include "progress_bar.h"

int main(int argc, char *argv[]){
	std::string output_basename;
	int num_bands_in_inverse = 8;
	int verbosity = 0;
	std::string mnf_statistics_basename;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"output-basename",		required_argument,	0,	'o'},
		{"mnf-statistics",		required_argument,	0,	'm'},
		{"verbose",			no_argument,		0,	'v'},
		{"num-bands-in-inverse",	required_argument,	0,	'b'},
		{0, 0, 0, 0}
	};
	char short_options[] = "hvo:m:b:";
	const char *option_descriptions[] = {
		"Show help",
		"Output filename for compressed image (mandatory)",
		"MNF statistics base filename (mandatory)",
		"Verbose output",
		"Number of bands in inverse, i.e. how compressed the image should be (default: 8)",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n"
			    "Will assume stdin input when no input filename is defined. \n";

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
			case 'm':
				mnf_statistics_basename = std::string(optarg);
			break;
			case 'v':
				verbosity++;
			break;
			case 'b': //number of bands in inverse
				num_bands_in_inverse = atoi(optarg);
			break;
		}
		if (flag == -1){
			break;
		}
	}

	//input filename
	std::string input_filename;
	if (optind < argc) {
		input_filename = std::string(argv[optind]);
	}

	//output filename
	if (output_basename.empty() || mnf_statistics_basename.empty()) {
		fprintf(stderr, "Output basename and mnf statistics input name must be defined.\n");
		exit(1);
	}

	//read hyperspectral header
	struct hyspex_header header;
	hyperspectral_err_t errcode;
	errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	//prepare hyperspectral image reader
	int stdin_options = STDIN_NO_ACCUMULATION;
	hyperspectral_reader_t *hyperspectral_reader = hyperspectral_create_reader(header, input_filename, stdin_options);

	//prepare MNF statistics
	mnf_err_t mnf_errcode;
	mnf_statistics_t *statistics = mnf_statistics_create(header);
	mnf_errcode = mnf_statistics_from_file(mnf_statistics_basename, header, statistics);
	if (mnf_errcode != MNF_NO_ERR) {
		fprintf(stderr, "Error reading statistics files with basename %s: %s\n", mnf_statistics_basename.c_str(), mnf_error_message(mnf_errcode));
		exit(1);
	}

	//compress image
	hyperspectral_compression_t *compressed_image = hyperspectral_create_compression(statistics, hyperspectral_reader, num_bands_in_inverse, verbosity);

	//write compression to file
	hyperspectral_compression_to_file(output_basename + ".hcompr", compressed_image);

	hyperspectral_free_compression(&compressed_image);
}

