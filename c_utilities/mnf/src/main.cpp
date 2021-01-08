#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "mnf.h"
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>
#include "getopt_long_helpers.h"
#include "progress_bar.h"
#include "string_helpers.h"

#define OPT_FORWARD_ONLY 205
#define OPT_INVERSE_ONLY 206
#define OPT_STATISTICS_FILE 207

int main(int argc, char *argv[]){
	int num_bands_in_inverse = 8;
	bool force_statistics_recalculation = false;
	int verbosity = 0;
	bool write_forward_mnf = true;
	bool write_inverse_mnf = true;
	bool output_to_stdout = true;

	std::string statistics_basename;
	std::string hyperspectral_basename;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"hyperspectral-file-output",	optional_argument,	0,	'w'},
		{"statistics-only",		no_argument,		0,	's'},
		{"statistics-file",		required_argument,	0,	OPT_STATISTICS_FILE},
		{"force-statistics-calculation",no_argument,		0,	'f'},
		{"num-bands-in-inverse",	required_argument,	0,	'b'},
		{"verbose",			no_argument,		0,	'v'},
		{"forward-only", 		no_argument,		0,	OPT_FORWARD_ONLY},
		{"inverse-only",		no_argument,		0,	OPT_INVERSE_ONLY},
		{0, 0, 0, 0}
	};
	char short_options[] = "";
	const char *option_descriptions[] = {
		"Show help",
		"Write hyperspectral output to file instead of stdout",
		"Write only statistics to file",
		"Specify statistics output basename",
		"Force recalculation of image and noise statistics even if the file already exists",
		"Number of bands to use in the inverse transform (default: 8)",
		"Increase verbosity",
		"Write only forward transform (disables inverse transform)",
		"Write only inverse transform (disables forward transform)",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n"
			    "Will assume stdin input when no input filename is defined. "
			    "Output is printed directly to stdout unless hyperspectral "
			    "file output is switched on.\n";

	//parse input arguments
	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;
			case 's': //generate only statistics
				write_forward_mnf = false;
				write_inverse_mnf = false;
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
			case 'w':
				output_to_stdout = false;
				if (optarg) {
					hyperspectral_basename = std::string(optarg);
				}
			break;
			case OPT_STATISTICS_FILE:
				statistics_basename = std::string(optarg);
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

	//define hyperspectral output filenames
	if (hyperspectral_basename.empty() && !input_filename.empty()) {
		hyperspectral_basename = get_basename(input_filename) + "_mnf";
	} else if (hyperspectral_basename.empty() && !output_to_stdout) {
		fprintf(stderr, "Input arrives from stdin and output will be written to file, but no output filename defined.\n");
		exit(1);
	}

	std::string mnf_forward_filename;
	std::string mnf_inverse_filename;
	if (!hyperspectral_basename.empty()) {
		mnf_forward_filename = hyperspectral_basename + "_forwardtransformed";
		mnf_inverse_filename = hyperspectral_basename + "_inversetransformed";
	}

	if (output_to_stdout && write_forward_mnf && write_inverse_mnf) {
		fprintf(stderr, "Cannot write both forward and inverse transform to stdout. Specify either of those transforms, or specify that output should be written to a file (see --help).\n");
		exit(1);
	}

	//define statistics output filename
	if (statistics_basename.empty() && !hyperspectral_basename.empty()) {
		statistics_basename = hyperspectral_basename;
	} else if (statistics_basename.empty()) {
		fprintf(stderr, "Need to define statistics output basename when using stdin as hyperspectral input (see --help).\n");
		exit(1);
	}

	//read hyperspectral header
	struct hyspex_header header;
	hyperspectral_err_t errcode;
	if (input_filename.empty()) {
		if (verbosity > 0) fprintf(stderr, "Reading image from stdin\n");
		errcode = hyperspectral_read_header_from_stdin(&header);
	} else {
		if (verbosity > 0) fprintf(stderr, "Reading image from file\n");
		errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	}
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	bool collect_statistics = !mnf_statistics_exist(statistics_basename) || force_statistics_recalculation;

	//prepare hyperspectral image reader
	int stdin_options = STDIN_NO_ACCUMULATION;
	if (input_filename.empty() && collect_statistics && (write_forward_mnf || write_inverse_mnf)) {
		fprintf(stderr, "Reading image from stdin and both collecting statistics and transforming the data, will accumulate the data in an array. Warning: Will consume memory equal to the total size of the hyperspectral input.\n");
		stdin_options = STDIN_ACCUMULATE_DATA;
	}
	hyperspectral_reader_t *hyperspectral_reader = hyperspectral_create_reader(header, input_filename, stdin_options);

	//prepare MNF statistics
	mnf_err_t mnf_errcode;
	mnf_statistics_t *statistics = mnf_statistics_create(header);
	mnf_transform_t *transform = mnf_transform_create(header, num_bands_in_inverse);
	float *image_line = new float[header.samples*header.bands]();
	struct hyspex_header line_header = header;
	line_header.lines = 1;

	if (!collect_statistics) {
		if (verbosity > 0) fprintf(stderr, "Reading MNF statistics from file\n");

		mnf_errcode = mnf_statistics_from_file(statistics_basename, header, statistics);
		if (mnf_errcode != MNF_NO_ERR) {
			fprintf(stderr, "Error reading statistics files with basename %s: %s\n", statistics_basename.c_str(), mnf_error_message(mnf_errcode));
			exit(1);
		}
	} else {
		for (int i=0; i < header.lines; i++) {
			if (verbosity > 0) print_progress_bar("Calculate MNF statistics", i, i-1, header.lines-1);

			errcode = hyperspectral_line_data_float(hyperspectral_reader, i, image_line);
			if (errcode != HYPERSPECTRAL_NO_ERR) {
				fprintf(stderr, "Image reading failed: %s\n", hyperspectral_error_message(errcode));
				exit(1);
			}
			mnf_update_statistics(statistics, line_header, image_line);
		}

		if (!statistics_basename.empty()) {
			if (verbosity > 0) fprintf(stderr, "Writing MNF statistics to file\n");
			mnf_statistics_to_file(statistics_basename, header, statistics);
		}
	}

	if (write_forward_mnf || write_inverse_mnf) {
		//prepare mnf transform
		mnf_errcode = mnf_transform_calculate(header, statistics, transform);
		if (mnf_errcode != MNF_NO_ERR) {
			fprintf(stderr, "Error calculating MNF transform: %s\n", mnf_error_message(mnf_errcode));
			exit(1);
		}

		//prepare output streams
		size_t data_size = header.bands*header.samples*header.lines*sizeof(float);
		hyperspectral_file_t mnf_forward_file;
		hyperspectral_file_t mnf_inverse_file;
		if (output_to_stdout) {
			mnf_forward_file = hyperspectral_open_stdout_file(header.lines, header.bands, header.samples, header.wlens);
			mnf_inverse_file = mnf_forward_file;
		} else {
			if (write_forward_mnf) {
				mnf_forward_file = hyperspectral_open_write_file(mnf_forward_filename, header.bands, header.samples, header.wlens);
			}
			if (write_inverse_mnf) {
				mnf_inverse_file = hyperspectral_open_write_file(mnf_inverse_filename, header.bands, header.samples, header.wlens);
			}
		}

		//run forward/inverse MNF and write to file(s)
		for (int i=0; i < header.lines; i++) {
			if (verbosity > 0) print_progress_bar("Applying MNF transform to data", i, i-1, header.lines-1);

			//acquire data
			size_t num_bytes = 0;
			size_t num_bytes_in_line = header.samples*header.bands*sizeof(float);
			errcode = hyperspectral_line_data_float(hyperspectral_reader, i, image_line);
			if (errcode != HYPERSPECTRAL_NO_ERR) {
				fprintf(stderr, "Image reading failed: %s\n", hyperspectral_error_message(errcode));
				exit(1);
			}

			//forward mnf
			if (write_forward_mnf) {
				mnf_run_forward(transform, line_header, image_line);
				num_bytes = hyperspectral_write_to_file(&mnf_forward_file, num_bytes_in_line, (char*)image_line);
				if (num_bytes != num_bytes_in_line) {
					fprintf(stderr, "Output stream broke off.\n");
					break;
				}
			}

			//inverse mnf
			if (write_inverse_mnf && !write_forward_mnf) {
				mnf_run_full(transform, line_header, image_line);
			} else if (write_inverse_mnf) {
				mnf_run_inverse(transform, line_header, image_line);
			}

			if (write_inverse_mnf) {
				num_bytes = hyperspectral_write_to_file(&mnf_inverse_file, num_bytes_in_line, (char*)image_line);
				if (num_bytes != num_bytes_in_line) {
					fprintf(stderr, "Output stream broke off.\n");
					break;
				}
			}
		}

		if (write_forward_mnf) hyperspectral_close_write_file(&mnf_forward_file);
		if (write_inverse_mnf) hyperspectral_close_write_file(&mnf_inverse_file);
	}

	mnf_free_transform(&transform);
	mnf_free_statistics(&statistics);
}
