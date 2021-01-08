#include <hyperspectral/readimage.h>
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
const int REQUIRED_NUM_COORDINATES = 2;

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
		"Specify coordinates from which to extract spectra, e.g. -c=LINE,SAMPLE",
		"Do not display wavelengths in the first printed column"
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [options] HYPERSPECTRAL_FILE\n";

	//discrete coordinates
	std::vector<std::string> input_coordinates;

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
				input_coordinates.push_back(std::string(optarg));
			break;
			case OPT_WITHOUT_WAVELENGTHS:
				display_wavelengths = false;
			break;
		}
	}

	//input sanity checks
	if (input_coordinates.size() == 0) {
		std::cerr << "No image coordinates specified, exiting." << std::endl;
		return -1;
	}

	//parse input coordinates
	std::vector<std::pair<int, int> > coordinates;
	for (int i=0; i < input_coordinates.size(); i++) {
		std::vector<int> input_coordinate = split_string(input_coordinates[i], ",");
		if (input_coordinate.size() != REQUIRED_NUM_COORDINATES) {
			std::cerr << "Error in coordinate " << input_coordinates[i] << ": Not the required number of coordinates (" << REQUIRED_NUM_COORDINATES << ")" << std::endl;
			exit(1);
		}

		coordinates.push_back(std::pair<int,int>(input_coordinate[0], input_coordinate[1]));
	}
	std::sort(coordinates.begin(), coordinates.end());

	//get input filename
	std::string input_filename;
	bool read_from_stdin = false;
	if (optind >= argc) {
		//no input files defined
		read_from_stdin = true;
		input_filename = "STDIN"; //for debug printing
	} else {
		input_filename = std::string(argv[optind]);
	}

	//read hyperspectral image
	struct hyspex_header header;
	hyperspectral_image_t image;
	hyperspectral_err_t errcode;

	if (!read_from_stdin) {
		errcode = hyperspectral_read_header_and_image(input_filename.c_str(), &header, &image);
	} else {
		errcode = hyperspectral_read_header_from_stdin(&header);
	}

	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error in reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	//collect spectra in line-order
	float *spectra = new float[header.bands*coordinates.size()]();
	char *read_line_data;
	if (read_from_stdin) {
		read_line_data = (char*)malloc(sizeof(float)*header.bands*header.samples);
	}
	int curr_stdin_line = -1;
	for (int i=0; i < coordinates.size(); i++) {
		int line = coordinates[i].first;
		int sample = coordinates[i].second;

		if (read_from_stdin) {
			//read from stdin until we reach the line we want
			while (line > curr_stdin_line) {
				int num_elements = fread(read_line_data, hyperspectral_element_size(header), header.bands*header.samples, stdin);
				if (num_elements < header.bands*header.samples) {
					fprintf(stderr, "Error in stdin reading.\n");
					exit(1);
				}
				curr_stdin_line++;
			}
		} else {
			read_line_data = hyperspectral_char_line(&image, line);
		}

		for (int j=0; j < header.bands; j++) {
			float value = 0;
			if (header.datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
				value = ((float*)read_line_data)[j*header.samples + sample];
			} else if (header.datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
				value = ((uint16_t*)read_line_data)[j*header.samples + sample];
			}
			spectra[j*coordinates.size() + i] = value;
		}
	}

	if (read_from_stdin) {
		fclose(stdin);
		free(read_line_data);
	}

	//output all spectra to stdout
	for (int j=0; j < header.bands; j++) {
		if (display_wavelengths) {
			std::cout << header.wlens[j] << " ";
		}
		for (int i=0; i < coordinates.size(); i++) {
			std::cout << spectra[j*coordinates.size() + i] << " ";
		}
		std::cout << std::endl;
	}

	hyperspectral_free_data(&image);
}
