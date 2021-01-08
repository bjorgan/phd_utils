#include <iostream>
#include <string>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "mnf.h"
#include <hyperspectral/readimage.h>
#include "getopt_long_helpers.h"

bool hyperspectral_read_line_from_stdin(struct hyspex_header header, char *line);

void hyperspectral_line_to_float(struct hyspex_header header, char *line, float *output);

int main(int argc, char *argv[]){
	bool output_to_stdout = true;

	bool collect_statistics = false;
	bool write_forward_transform = false;
	bool write_inverse_transform = false;
	bool write_statistics = false;

	int num_bands_in_inverse = 0;

	std::string statistics_basename;
	std::string hyperspectral_basename;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"output-to-file",		required_argument,	0,	'o'},
		{0, 0, 0, 0}
	};
	char short_options[] = "h:o";
	const char *option_descriptions[] = {
		"Show help",
		"Output statistics and image to file",
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
		}
		if (flag == -1){
			break;
		}
	}

	//input filename/stdin
	std::string input_filename;
	bool read_from_stdin = false;
	if (optind >= argc) {
		input_filename = "STDIN"; //for debug printing
		read_from_stdin = true;
	} else {
		input_filename = std::string(argv[optind]);
	}

	//output filename
	/*
	if (output_basename.empty() && !input_filename.empty()) {
		char *input_copy = strdup(input_filename.c_str());
		char *basename_str = basename(input_copy);
		output_basename = std::string(basename_str);
		output_basename += "_mnf";
		free(input_copy);
	}
	*/

	std::string mnf_forward_filename = hyperspectral_basename + "_forwardtransformed";
	std::string mnf_inverse_filename = hyperspectral_basename + "_inversetransformed";

	if (!collect_statistics) {
		write_statistics = false;
	}

	if (write_statistics && statistics_basename.empty()) {
		fprintf(stderr, "Need output name for writing statistics.\n");
		exit(1);
	}

	if (!collect_statistics && statistics_basename.empty()) {
		fprintf(stderr, "Need input name for reading statistics.\n");
		exit(1);
	}

	if (read_from_stdin && collect_statistics && (write_forward_transform || write_inverse_transform)) {
		fprintf(stderr, "Need to run two passes on the data, but reading from stdin.\n");
		exit(1);
	}

	if (collect_statistics && (write_forward_transform || write_inverse_transform)) {
		fprintf(stderr, "Need to run two passes on the data, but not implemented here.\n");
		exit(1);
	}

	if (output_to_stdout && write_forward_transform && write_inverse_transform) {
		fprintf(stderr, "Cannot write both forward and inverse transform to stdout.");
		exit(1);
	}

	//read hyperspectral image
	struct hyspex_header header;
	hyperspectral_image_t hyperspectral_image;
	hyperspectral_err_t errcode;
	if (read_from_stdin) {
		errcode = hyperspectral_read_header_from_stdin(&header);
	} else {
		errcode = hyperspectral_read_header_and_image(input_filename.c_str(), &header, &hyperspectral_image);
	}
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	float *float_line_data = new float[header.bands*header.samples]();
	char *char_line_data = NULL;
	if (read_from_stdin) {
		char_line_data = new char[hyperspectral_element_size(header)*header.samples*header.bands];
	}

	//MNF transform
	mnf_statistics_t statistics = mnf_statistics_create(header);
	mnf_transform_t transform = mnf_transform_create(header, num_bands_in_inverse);
	struct hyspex_header line_header = header;
	line_header.lines = 1;

	//prepare hyperspectral output streams
	hyperspectral_file_t mnf_forward_file, mnf_inverse_file;
	if (output_to_stdout && (write_forward_transform || write_inverse_transform)) {
		mnf_forward_file = hyperspectral_open_stdout_file(header.lines, header.bands, header.samples, header.wlens);
		mnf_inverse_file = mnf_forward_file;
	} else if (!output_to_stdout) {
		if (write_forward_transform) {
			mnf_forward_file = hyperspectral_open_write_file(mnf_forward_filename, header.bands, header.samples, header.wlens);
		}

		if (write_inverse_transform) {
			mnf_inverse_file = hyperspectral_open_write_file(mnf_inverse_filename, header.bands, header.samples, header.wlens);
		}
	}

	if (collect_statistics) {
		for (int line_number=0; line_number < header.lines; line_number++) {
			//prepare float array
			if (read_from_stdin) {
				if (!hyperspectral_read_line_from_stdin(header, char_line_data)) {
					fprintf(stderr, "stdin stream broke off unexpectedly.\n");
					exit(1);
				}
			} else {
				char_line_data = hyperspectral_char_line(&hyperspectral_image, line_number);
			}
			hyperspectral_line_to_float(header, char_line_data, float_line_data);
			size_t num_bytes_in_line = header.samples*header.bands*sizeof(float);

			//update statistics
			mnf_update_statistics(&statistics, line_header, float_line_data);
		}
	}

	if (write_forward_transform || write_inverse_transform) {
		for (int line_number=0; line_number < header.lines; line_number++) {
			//prepare float array
			if (read_from_stdin) {
				if (!hyperspectral_read_line_from_stdin(header, char_line_data)) {
					fprintf(stderr, "stdin stream broke off unexpectedly.\n");
					exit(1);
				}
			} else {
				char_line_data = hyperspectral_char_line(&hyperspectral_image, line_number);
			}
			hyperspectral_line_to_float(header, char_line_data, float_line_data);
			size_t num_bytes_in_line = header.samples*header.bands*sizeof(float);

			//run forward/inverse transforms and write to files appropriately
			if (write_forward_transform && !write_inverse_transform) {
				mnf_run_forward(&transform, line_header, float_line_data);
				hyperspectral_write_to_file(&mnf_forward_file, num_bytes_in_line, (char*)float_line_data);
			} else if (write_inverse_transform && !write_forward_transform) {
				mnf_run_full(&transform, line_header, float_line_data);
				hyperspectral_write_to_file(&mnf_inverse_file, num_bytes_in_line, (char*)float_line_data);
			} else {
				mnf_run_forward(&transform, line_header, float_line_data);
				hyperspectral_write_to_file(&mnf_forward_file, num_bytes_in_line, (char*)float_line_data);
				mnf_run_inverse(&transform, line_header, float_line_data);
				hyperspectral_write_to_file(&mnf_inverse_file, num_bytes_in_line, (char*)float_line_data);
			}
		}
	}

	//close files
	if (write_forward_transform) {
		hyperspectral_close_write_file(&mnf_forward_file);
	}

	if (write_inverse_transform) {
		hyperspectral_close_write_file(&mnf_inverse_file);
	}
}

bool hyperspectral_read_line_from_stdin(struct hyspex_header header, char *line)
{
	int num_elements = fread(line, hyperspectral_element_size(header), header.samples*header.bands, stdin);
	return num_elements == header.samples*header.bands;
}

void hyperspectral_line_to_float(struct hyspex_header header, char *line, float *output)
{
	if (header.datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
		memcpy(output, line, sizeof(float)*header.samples*header.bands);
	} else if (header.datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
		for (int i=0; i < header.samples*header.bands; i++) {
			output[i] = ((uint16_t*)line)[i];
		}
	}
}
