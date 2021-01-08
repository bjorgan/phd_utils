#include <hyperspectral/readimage_stdin.h>
#include <iostream>
#include <getopt.h>
#include <string>
#include <cstring>
#include <libgen.h>
#include <unistd.h>
#include "getopt_long_helpers.h"
#include <algorithm>
#include "string_helpers.h"

#define OPT_WITHOUT_WAVELENGTHS 203

///Required number of coordinates to specify on command line in format y,x (2)
const int REQUIRED_NUM_COORDINATES = 4;

size_t hyperspectral_element_size(struct hyspex_header header); //from readimage.cpp

int main(int argc, char *argv[])
{
	//command line options
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"coordinate",			required_argument,	0,	'c'},
		{"no-wavelengths",		no_argument,		0,	OPT_WITHOUT_WAVELENGTHS},
		{0, 0, 0, 0}
	};
	char short_options[] = "c:h";
	const char* option_descriptions[] = {
		"Show help",
		"Specify coordinate from which to extract a mean spectrum, e.g. -c=STARTLINE:ENDLINE,STARTSAMPLE:ENDSAMPLE",
		"Do not display wavelengths in the first printed column"
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [options] HYPERSPECTRAL_FILE\n";

	//discrete coordinates
	std::string input_coordinates;

	//whether wavelengths should be printed to stdout
	bool display_wavelengths = true;

	while (1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, short_options, long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'h':
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				return 0;
			break;
			case 'c':
				input_coordinates = std::string(optarg);
			break;
			case OPT_WITHOUT_WAVELENGTHS:
				display_wavelengths = false;
			break;
		}
	}

	//input sanity checks
	if (input_coordinates.empty()) {
		std::cerr << "No image coordinates specified, exiting." << std::endl;
		return -1;
	}

	//parse input coordinates
	std::vector<int> coordinates = split_string(input_coordinates, ",:");
	if (coordinates.size() != REQUIRED_NUM_COORDINATES) {
		std::cerr << "Error in coordinate " << input_coordinates << ": Not the required number of coordinates (" << REQUIRED_NUM_COORDINATES << ")" << std::endl;
		exit(1);
	}
	int start_line = coordinates[0];
	int end_line = coordinates[1];
	int start_sample = coordinates[2];
	int end_sample = coordinates[3];

	//get input filename
	std::string input_filename;
	bool read_from_stdin = false;
	if (optind < argc) {
		input_filename = std::string(argv[optind]);
	}

	//read hyperspectral image
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error in reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	hyperspectral_reader_t *reader = hyperspectral_create_reader(header, input_filename);

	//acquire mean spectrum
	int num_spectra = 1;
	double *spectra = new double[header.bands*num_spectra]();
	float *line = new float[header.bands*header.samples]();
	for (int i=start_line; i < end_line; i++) {
		errcode = hyperspectral_line_data_float(reader, i, line);
		if (errcode) {
			fprintf(stderr, "Error in reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
			exit(1);
		}

		for (int k=0; k < header.bands; k++) {
			for (int j=start_sample; j < end_sample; j++) {
				spectra[k] += line[k*header.samples + j];
			}
		}
	}
	for (int i=0; i < header.bands; i++) {
		spectra[i] /= 1.0*(end_line - start_line)*(end_sample - start_sample);
	}

	//output all spectra to stdout
	for (int j=0; j < header.bands; j++) {
		if (display_wavelengths) {
			std::cout << header.wlens[j] << " ";
		}
		for (int i=0; i < num_spectra; i++) {
			std::cout << spectra[j*num_spectra + i] << " ";
		}
		std::cout << std::endl;
	}
}
