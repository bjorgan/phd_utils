//==============================================================================
// Copyright 2017 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>
#include "getopt_long_helpers.h"
#include "string_helpers.h"

#define OPT_START_LINE 205
#define OPT_END_LINE 206
#define OPT_START_SAMPLE 207
#define OPT_END_SAMPLE 208
#define OPT_START_BAND 209
#define OPT_END_BAND 210
#define OPT_SAMPLES 211
#define OPT_LINES 212
#define OPT_BANDS 213

/**
 * Rectify ranges of image subset against image header.
 *
 * \param header Image header
 * \param subset Subset
 * \return Rectified image subset
 **/
struct image_subset rectify_subset(struct hyspex_header header, struct image_subset subset);

void get_start_end(std::string input, int *ret_start, int *ret_end)
{
	std::vector<int> coordinates = split_string(input, ":");
	if (coordinates.size() != 2) {
		fprintf(stderr, "Insufficient number of arguments: %s\n", input.c_str());
		exit(1);
	}
	*ret_start = coordinates[0];
	*ret_end = coordinates[1];
}

int main(int argc, char *argv[]){
	bool output_to_stdout = true;

	std::string hyperspectral_basename;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"hyperspectral-file-output",	optional_argument,	0,	'w'},
		{"startline",			required_argument,	0,	OPT_START_LINE},
		{"endline",			required_argument,	0,	OPT_END_LINE},
		{"startsample",			required_argument,	0,	OPT_START_SAMPLE},
		{"endsample",			required_argument,	0,	OPT_END_SAMPLE},
		{"startband",			required_argument,	0,	OPT_START_BAND},
		{"endband",			required_argument,	0,	OPT_END_BAND},
		{"samples",			required_argument,	0,	OPT_SAMPLES},
		{"lines",			required_argument,	0,	OPT_LINES},
		{"bands",			required_argument,	0,	OPT_BANDS},
		{0, 0, 0, 0}
	};
	char short_options[] = "";
	const char *option_descriptions[] = {
		"Show help",
		"Write hyperspectral output to file instead of stdout",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n"
			    "Will assume stdin input when no input filename is defined. "
			    "Output is printed directly to stdout unless hyperspectral "
			    "file output is switched on.\n";

	struct image_subset subset = {0};

	bool separate_coordinates_specified = false;
	bool coordinate_ranges_specified = false;

	//parse input arguments
	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;
			case 'w':
				output_to_stdout = false;
				if (optarg) {
					hyperspectral_basename = std::string(optarg);
				}
			break;

			case OPT_START_LINE:
				subset.start_line = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_END_LINE:
				subset.end_line = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_START_SAMPLE:
				subset.start_sample = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_END_SAMPLE:
				subset.end_sample = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_START_BAND:
				subset.start_band = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_END_BAND:
				subset.end_band = atoi(optarg);
				separate_coordinates_specified = true;
			break;

			case OPT_SAMPLES:
				get_start_end(optarg, &subset.start_sample, &subset.end_sample);
				coordinate_ranges_specified = true;
			break;

			case OPT_LINES:
				get_start_end(optarg, &subset.start_line, &subset.end_line);
				coordinate_ranges_specified = true;
			break;

			case OPT_BANDS:
				get_start_end(optarg, &subset.start_band, &subset.end_band);
				coordinate_ranges_specified = true;
			break;

		}
		if (flag == -1){
			break;
		}
	}

	if (coordinate_ranges_specified && separate_coordinates_specified) {
		fprintf(stderr, "Image range specified using both separate specifier and a combined specifier, this is incompatible.\n");
		exit(1);
	}

	//input filename
	std::string input_filename;
	if (optind < argc) {
		input_filename = std::string(argv[optind]);
	}

	//define hyperspectral output filenames
	if (hyperspectral_basename.empty() && !input_filename.empty()) {
		hyperspectral_basename = get_basename(input_filename) + "_cropped";
	} else if (hyperspectral_basename.empty() && !output_to_stdout) {
		fprintf(stderr, "Input arrives from stdin and output will be written to file, but no output filename defined.\n");
		exit(1);
	}
	std::string output_filename;
	if (!output_to_stdout) {
		output_filename = hyperspectral_basename;
	}

	//read hyperspectral header
	struct hyspex_header header;
	hyperspectral_err_t errcode;
	errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	subset = rectify_subset(header, subset);

	//prepare hyperspectral image reader
	hyperspectral_reader_t *reader = hyperspectral_create_reader(header, input_filename);

	int start_line = subset.start_line;

	struct hyspex_header output_header = hyperspectral_header_from_subset(header, subset);
	output_header.offset = 0;
	size_t element_size = hyperspectral_element_size(header);
	char *output_line = new char[element_size*output_header.lines*output_header.bands*output_header.samples]();

	hyperspectral_file_t writer = hyperspectral_open(output_filename, output_header);

	//crop image
	for (int i=start_line; i < subset.end_line; i++) {
		const char *line = hyperspectral_line_data(reader, i);
		int k = 0;
		for (int j=subset.start_band; j < subset.end_band; j++) {
			memcpy(output_line + k*output_header.samples*element_size, line + j*header.samples*element_size + subset.start_sample*element_size, element_size*output_header.samples);
			k++;
		}
		size_t line_size = element_size*output_header.samples*output_header.bands;
		size_t written_size = hyperspectral_write_to_file(&writer, line_size, output_line);
		if (line_size != written_size) {
			fprintf(stderr, "Output stream broke off.\n");
			break;
		}
	}
	hyperspectral_close_write_file(&writer);
}

struct image_subset rectify_subset(struct hyspex_header header, struct image_subset subset)
{
	if (!subset.end_sample || (subset.end_sample > header.samples)) subset.end_sample = header.samples;
	if (!subset.end_line || (subset.end_line > header.lines)) subset.end_line = header.lines;
	if (!subset.end_band || (subset.end_band > header.bands)) subset.end_band = header.bands;

	if (subset.start_sample >= subset.end_sample) subset.start_sample = subset.end_sample-1;
	if (subset.start_line >= subset.end_line) subset.start_line = subset.end_line-1;
	if (subset.start_band >= subset.end_band) subset.start_band = subset.end_band-1;

	return subset;
}
