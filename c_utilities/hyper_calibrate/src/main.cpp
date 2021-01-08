#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include "progress_bar.h"
#include "getopt_long_helpers.h"
#include "string_helpers.h"
#include "calibration.h"
#include "calibration_io.h"
#include <unistd.h>
#include <hyperspectral/readimage_stdin.h>
#include "divide_by_max.h"
using namespace std;

#define OPT_NO_RGB_IMAGE 203
#define OPT_NO_HYPERSPECTRAL 204
#define OPT_DIVIDE_BY_MAX 205
#define OPT_CROP_IN_SAMPLE_DIRECTION 206
#define OPT_INTEGRATE_ALL 207
#define OPT_REFL_STD_LINES 208
#define OPT_REFL_STD_SAMPLES 209
#define OPT_REFL_STD_CONSTANT_SPEC 210

std::string remove_file_ending(std::string filename); //from readimage

int main(int argc, char *argv[])
{
	//search direction for locating the reflectance standard
	search_direction direction = START_FROM_START;

	//output filename (assumed to be explicitly set if non-empty)
	string out_filename;

	//filename for reflectance standard specification file (assumed to be explicitly set if non-empty)
	string specs_filename;

	//print verbosity
	int verbosity = 0;

	//specified end of reflectance standard (assumed to not be set if non-zero. Will be offset from start or end of image depending on search direction.)
	int specified_cal_end = 0;

	//find and integrate reflectance standard
	bool integrate_reflectance_standard = true;

	//alternative path to calibration info if above is set to false
	std::string calibration_info_file;

	//generate RGB image and indicate position of reflectance standard
	bool generate_rgb_image = true;

	//output to disk
	bool output_to_file = false;

	//output to standard output
	bool output_to_stdout = true;

	//do only reflectance standard integration
	bool integrate_only = false;

	//force reintegration of reflectance standard even if calibration info file already exists
	bool force_reintegration = false;

	//calibrate by the maximum value
	bool calibrate_by_maximum_value = false;

	//crop reflectance standard in sample direction using same algorithm as for the line direction cropping
	bool crop_in_sample_direction = false;

	//whether we should integrate over the entire reflectance standard or integrate sample by sample and keep spatial distribution
	bool integrate_all = false;

	//for optionally explicitly specifying reflectance standard positions if we don't want to search for it.
	struct image_subset refl_standard_positions_from_argv = {-1, -1, -1, -1};
	bool refl_standard_positions_specified_from_argv = false;

	//constant standard adjustment values, as an alternative to specification files
	bool adjust_standard_with_constant_value = false;
	double standard_adjustment_value = 1.0;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"output-to-file",		optional_argument,	0,	'o'},
		{"specs-filename",		required_argument,	0,	's'},
		{"search-from-end",		no_argument,		0,	'e'},
		{"cal-end-line",		required_argument,	0,	'l'},
		{"verbose",			no_argument,		0,	'v'},
		{"calibration-file",		required_argument,	0,	'c'},
		{"no-rgb",			no_argument,		0,	OPT_NO_RGB_IMAGE},
		{"force-reintegration",		no_argument,		0,	'f'},
		{"integrate-only",		no_argument,		0,	'i'},
		{"divide-by-max",		no_argument,		0,	OPT_DIVIDE_BY_MAX},
		{"crop-in-sample-direction",	no_argument,		0,	OPT_CROP_IN_SAMPLE_DIRECTION},
		{"integrate-all",		no_argument,		0,	OPT_INTEGRATE_ALL},
		{"force-refl-std-lines",	required_argument,	0,	OPT_REFL_STD_LINES},
		{"force-refl-std-samples",	required_argument,	0,	OPT_REFL_STD_SAMPLES},
		{"adjust-standard-with-constant-value",	required_argument,	0,	OPT_REFL_STD_CONSTANT_SPEC},
		{0, 0, 0, 0}
	};
	char short_options[] = "ho::s:el:vc:fi";
	const char *option_descriptions[] = {
		"Show help",
		"Output hyperspectral image to file instead of stdout. Can specify base filename for output, will otherwise assume working directory as output directory, with ${input_filename}_calibration.img and .hdr, and _calibration.calibration_info for metainfo",
		"Filename of specification file for what the reflectance from the reflectance standard is. It is assumed that the specification file consists of two columns, wavelengths and values between 0 and 1",
		"Search for reflectance standard from the end of the image",
		"Specify end line of reflectance standard (if searched from start), or offset counted from end (if searched from end)",
		"Increase verbosity of output",
		"Specify calibration info file and don't look for reflectance standard in image",
		"Don't output RGB image",
		"Reintegrate reflectance standard even if calibration info file already exists",
		"Do only reflectance standard integration and write metafile, do not calibrate hyperspectral image",
		"Divide by the maximum value along each spectrum instead of applying reflectance standard-based calibration. Ignores the rest of the options",
		"Try to crop along sample direction after having found reflectance standard in line directions",
		"Integrate over entire reflectance standard instead of retaining spatial distribution of the light source",
		"Instead of detecting reflectance standard positions automatically: Set start and end line explicitly. Assumes form startline:endline",
		"Same, but in sample direction. Specify both to set a rectangle",
		"Adjust reflectance standard with a constant value instead of a spectrum from a specification file",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n\nSee --force-refl-std* options for explicitly specifying reflectance standard positions. Otherwise defaults to searching for it from the start of the image, assuming FOV-covering reflectance standard. See other options below for perturbations to this behavior.\n";

	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);

		if (flag == -1){
			break;
		}

		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;

			case 'o': //output to file
				output_to_file = true;
				output_to_stdout = false;
				if (optarg) {
					out_filename = string(optarg);
				}
			break;

			case 's': //reflectance standard specification
				specs_filename = string(optarg);
			break;

			case 'e': //search direction
				direction = START_FROM_END;
			break;

			case 'l': //end of reflectance standard
				specified_cal_end = strtod(optarg, NULL);
			break;

			case 'v': //verbosity level
				verbosity++;
			break;

			case 'c': //calibration file
				calibration_info_file = std::string(optarg);
				integrate_reflectance_standard = false;
			break;

			case OPT_NO_RGB_IMAGE: //should output no RGB image
				generate_rgb_image = false;
			break;

			case 'f': //force reintegration even if file exists
				force_reintegration = true;
			break;

			case 'i': //integrate only
				integrate_only = true;
			break;

			case OPT_DIVIDE_BY_MAX:
				calibrate_by_maximum_value = true;
			break;

			case OPT_CROP_IN_SAMPLE_DIRECTION:
				crop_in_sample_direction = true;
			break;

			case OPT_INTEGRATE_ALL:
				integrate_all = true;
			break;

			case OPT_REFL_STD_LINES:
			{
				std::vector<int> lines = split_string(std::string(optarg), ":");
				refl_standard_positions_from_argv.start_line = lines[0];
				refl_standard_positions_from_argv.end_line = lines[1];
				refl_standard_positions_specified_from_argv = true;
			}
			break;
			case OPT_REFL_STD_SAMPLES:
			{
				std::vector<int> samples = split_string(std::string(optarg), ":");
				refl_standard_positions_from_argv.start_sample = samples[0];
				refl_standard_positions_from_argv.end_sample = samples[1];
				refl_standard_positions_specified_from_argv = true;
			}
			break;
			case OPT_REFL_STD_CONSTANT_SPEC:
				adjust_standard_with_constant_value = true;
				standard_adjustment_value = strtod(optarg, NULL);
				break;

			default:
				exit(1);
			break;
		}
	}

	if (integrate_only) {
		output_to_file = false;
		output_to_stdout = false;
	}

	//input filename
	string filename;
	if (optind < argc) {
		filename = string(argv[optind]);
	}

	if (output_to_file && out_filename.empty() && filename.empty()) {
		fprintf(stderr, "Need to define output filename.\n");
		exit(1);
	}

	//output filename
	if (out_filename.empty() && output_to_file){
		//use working directory as output directory
		out_filename = get_basename(filename) + "_calibrated";
	}

	if (!out_filename.empty() && !filename.empty() && (remove_file_ending(filename) == out_filename)) {
		fprintf(stderr, "Output filename the same as input filename, exiting.\n");
		exit(1);
	}

	if (calibrate_by_maximum_value) {
		//skip the rest of the calibration main, and read the file using the general stdin/file interface
		calibrate_by_max_spectral_value(filename, out_filename);
		return 0;
	}

	if (calibration_info_file.empty() && filename.empty()) {
		fprintf(stderr, "Need calibration info filename when running in stdin input mode.\n");
		exit(1);
	}

	if (calibration_info_file.empty()) {
		calibration_info_file = out_filename + CALIBRATION_INFO_FILE_POSTFIX;
	}

	//check if calibration info file already exists
	if ((access(calibration_info_file.c_str(), F_OK) == 0) && !force_reintegration) {
		integrate_reflectance_standard = false;
	}

	if (integrate_reflectance_standard && filename.empty()) {
		fprintf(stderr, "Cannot integrate reflectance standard in stdin input mode.\n");
		exit(1);
	}

	//read header
	struct hyspex_header header;
	hyperspectral_err_t ret = hyperspectral_read_header(filename.c_str(), &header);
	if (ret != HYPERSPECTRAL_NO_ERR){
		fprintf(stderr, "Error in reading image: %s, %s\n", filename.c_str(), hyperspectral_error_message(ret));
		exit(1);
	}


	hyperspectral_calibration_info_t calibration_info = {0};
	if (integrate_reflectance_standard) {
		//read image for reflectance standard integration
		hyperspectral_image_t hyperspectral_data;
		ret = hyperspectral_read_image(filename.c_str(), header, &hyperspectral_data);
		if (ret != HYPERSPECTRAL_NO_ERR){
			fprintf(stderr, "Error in reading image: %s, %s\n", filename.c_str(), hyperspectral_error_message(ret));
			exit(1);
		}
		const char *image_data = hyperspectral_char(&hyperspectral_data);

		//get position for reflectance standard
		struct image_subset refl_std_pos;
		if (refl_standard_positions_specified_from_argv) {
			refl_std_pos = refl_standard_positions_from_argv;
			if (refl_std_pos.start_line == -1) refl_std_pos.start_line = 0;
			if (refl_std_pos.end_line == -1) refl_std_pos.end_line = header.lines-1;
			if (refl_std_pos.start_sample == -1) refl_std_pos.start_sample = 0;
			if (refl_std_pos.end_sample == -1) refl_std_pos.end_sample = header.samples-1;
		} else {
			if (specified_cal_end == 0) {
				if (verbosity > 0) fprintf(stderr, "Finding positions of reflectance standard.\n");
				refl_std_pos = calibration_find_FOV_refl_standard(direction, header, image_data, verbosity);

				if (crop_in_sample_direction) {
					if (verbosity > 0) fprintf(stderr, "Cropping reflectance standard in across-track direction.\n");
					calibration_crop_refl_standard_in_sample_direction(refl_std_pos, header, image_data, verbosity);
				}

			} else {
				refl_std_pos.start_sample = 0;
				refl_std_pos.end_sample = header.samples;

				if (direction == START_FROM_START) {
					refl_std_pos.start_line = 0;
					refl_std_pos.end_line = specified_cal_end;
				} else {
					refl_std_pos.start_line = header.lines - specified_cal_end;
					refl_std_pos.end_line = header.lines - 1;
				}
			}
		}

		//get calibration array
		if (!integrate_all) {
			//integrating only in the line direction, keeping the sample distribution in order to rectify for light source variations across field of view. This is not
			//desirable unless the reflectance standard covers the entire FOV, display warning if positions seem to indicate that it doesn't.
			if ((refl_std_pos.start_line != 0) && (refl_std_pos.end_line != header.lines-1)) fprintf(stderr, "Warning: Reflectance standard not specified to cover entire FOV, but keeping spatial distribution over the reflectance standard and extrapolating the edges regardless. Consider to use --integrate-all.\n");
			calibration_info = calibration_integrate(header, image_data, refl_std_pos, verbosity);
		} else {
			//integrate over the entire reflectance standard
			calibration_info = calibration_integrate_all(header, image_data, refl_std_pos, verbosity);
		}

		//adjust against specifications
		if (!specs_filename.empty()){
			if (verbosity > 0) fprintf(stderr, "Adjusting calibration array against specification file.\n");
			calibration_adjust(specs_filename, header, &calibration_info);
		} else if (adjust_standard_with_constant_value) {
			if (verbosity > 0) fprintf(stderr, "Adjusting calibration array against specified constant value.\n");
			calibration_adjust_with_constant_value(standard_adjustment_value, header, &calibration_info);
		} else {
			fprintf(stderr, "Warning: Calibration standard specification file not specified. Reflectance values will not be strictly correct since spectral variation of the reflectance standard is not taken into account.\n");
		}

		//write calibration info to file
		hyperspectral_calibration_info_to_file(&calibration_info, calibration_info_file);

		//write RGB image
		#ifdef WITH_OPENCV
		if (generate_rgb_image) {
			float *rgb_data = new float[header.samples*header.lines*header.default_bands.size()]();
			float *float_line = new float[header.samples*header.bands]();
			for (int i=0; i < header.lines; i++) {
				if (verbosity > 0) print_progress_bar("Generate RGB image.", i-1, i, header.lines-1);
				const char *line_data = hyperspectral_char_line(&hyperspectral_data, i);
				hyperspectral_line_to_float(header, line_data, float_line);
				for (int k=0; k < header.default_bands.size(); k++) {
					float *rgb_band = rgb_data + k*header.samples*header.lines;
					int band = header.default_bands[k];
					for (int j=0; j < header.samples; j++) {
						int cal_index = band*header.samples + j;
						if (calibration_info.num_samples == 1) {
							cal_index = band;
						}

						rgb_band[i*header.samples + j] = float_line[band*header.samples + j]/calibration_info.calibration_array[cal_index];
					}
				}
			}
			delete [] float_line;

			//get and write rgb image
			cv::Mat rgb_image = calibration_generate_rgb_image(header, rgb_data, refl_std_pos);
			cv::imwrite(calibration_info_file + "_refl_position.png", rgb_image);
			delete [] rgb_data;
		}
		#endif

		hyperspectral_free_data(&hyperspectral_data);
	} else {
		if (verbosity > 0) fprintf(stderr, "Reading calibration info from file.\n");
		calibration_info = hyperspectral_calibration_info_from_file(calibration_info_file);
		if (!hyperspectral_calibration_info_valid(header, &calibration_info)) {
			fprintf(stderr, "Calibration info %s not valid, exiting.\n", calibration_info_file.c_str());
			exit(1);
		}
	}

	if (output_to_file || output_to_stdout) {
		//reopen hyperspectral image as a stream (so that we can accept either stdin or file input depending on whether filename is defined)
		hyperspectral_reader_t *reader = hyperspectral_create_reader(header, filename);

		//prepare output
		float *output_data = new float[header.bands*header.samples]();
		hyperspectral_file_t hyperspectral_file;
		if (output_to_file) {
			hyperspectral_file = hyperspectral_open_write_file(out_filename, header.bands, header.samples, header.wlens);
		} else {
			hyperspectral_file = hyperspectral_open_stdout_file(header.lines, header.bands, header.samples, header.wlens);
		}

		//calibrate image
		for (int i=0; i < header.lines; i++){
			const char *line_data = hyperspectral_line_data(reader, i);
			if (line_data == NULL) {
				fprintf(stderr, "Stream closed.\n");
				exit(1);
			}
			calibration_run_single_line(&calibration_info, header, line_data, output_data);

			if (verbosity > 0) print_progress_bar("Calibrate image and write results to file.", i-1, i, header.lines-1);

			//write image to disk
			int num_bytes = header.samples*header.bands*sizeof(float);
			int written_bytes = hyperspectral_write_to_file(&hyperspectral_file, num_bytes, (char*)output_data);
			if (written_bytes < num_bytes) {
				fprintf(stderr, "Stream closed.\n");
				exit(1);
			}
		}
		hyperspectral_close_write_file(&hyperspectral_file);
		delete [] output_data;
	}

	hyperspectral_calibration_info_free(&calibration_info);
}

