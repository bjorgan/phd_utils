#include <fstream>
#include "masking.h"
#include <iostream>
#include <sys/time.h>
#include "getopt_long_helpers.h"
#include "string_helpers.h"
#include "progress_bar.h"
#include <hyperspectral/readimage_stdin.h>
using namespace std;

#define OPT_OUTPUT_MASK 204
#define OPT_OUTPUT_HYPERSPECTRAL_FILE 205
#define OPT_OUTPUT_MASK_ONLY 206
#define OPT_SAM_THRESH 207

void check_reading_errors(std::string filename, hyperspectral_err_t errcode)
{
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}
}

int main(int argc, char *argv[])
{
	std::string output_basename;
	std::string output_hyperspectral_basename;
	std::string spectral_directory;
	std::vector<std::string> spectral_files;
	std::vector<std::pair<int, int> > spectral_coordinates;
	double sam_thresh = SAM_THRESH_DEFAULT;
	int verbosity = 0;
	bool write_hyperspectral_output = true;
	bool write_hyperspectral_file_to_stdout = true;
	bool write_mask_to_file = false;

	//getopt properties
	struct option long_options[] = {
		{"help",				no_argument,		0,	'h'},
		{"spectral-directory",			required_argument,	0,	'd'},
		{"spectral-file",			required_argument,	0,	'f'},
		{"spectral-coordinate",			required_argument,	0,	'c'},
		{"verbose",				no_argument,		0,	'v'},
		{"write-mask",				no_argument,		0,	OPT_OUTPUT_MASK},
		{"write-mask-only",			no_argument,		0,	OPT_OUTPUT_MASK_ONLY},
		{"write-hyperspectral-image-to-file",	no_argument,		0,	OPT_OUTPUT_HYPERSPECTRAL_FILE},
		{"output-basename",			required_argument,	0,	'o'},
		{"sam-thresh",				required_argument,	0,	OPT_SAM_THRESH},
		{0, 0, 0, 0}
	};
	char short_options[] = "h:d:f:c:v";
	const char *option_descriptions[] = {
		"Show help",
		"Specify directory for masking spectra",
		"Specify file for masking spectra. Can specify multiple files. Takes precedence over directories",
		"Specify coordinates (-c LINE,SAMPLE) for masking spectra. Can specify multiple coordinates. Takes precedence over files and directories",
		"Verbose output",
		"Write mask to plain text file",
		"Write only mask, no hyperspectral output",
		"Write hyperspectral image to file instead of stdout",
		"Specify base filename for output files",
		"Specify SAM threshold",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n";
	int index;
	std::vector<int> coordinates;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;

			case 'd':
				spectral_directory = std::string(optarg);
			break;

			case 'f':
				spectral_files.push_back(std::string(optarg));
			break;

			case 'c':
				coordinates = split_string(std::string(optarg), ",");
				spectral_coordinates.push_back(std::pair<int,int>(coordinates[0], coordinates[1]));
			break;

			case 'v':
				verbosity++;
			break;

			case OPT_OUTPUT_MASK_ONLY:
				write_hyperspectral_output = false;
			case OPT_OUTPUT_MASK:
				write_mask_to_file = true;
			break;

			case OPT_OUTPUT_HYPERSPECTRAL_FILE:
				write_hyperspectral_file_to_stdout = false;
			break;

			case 'o':
				output_basename = std::string(optarg);
			break;

			case OPT_SAM_THRESH:
				sam_thresh = strtod(optarg, NULL);
			break;
		}
		if (flag == -1){
			break;
		}
	}

	string filename;
	if (optind < argc) {
		filename = std::string(argv[optind]);
	}

	if (!filename.empty() && output_basename.empty()) {
		output_basename = get_basename(filename) + "_masking";
	}

	if ((write_mask_to_file || !write_hyperspectral_file_to_stdout) && output_basename.empty()) {
		fprintf(stderr, "Flags for writing mask or hyperspectral image to file are specified, but no output file, and input is arriving from stdin, so cannot generate output filename from input filename. Specify output basename using --output-basename.\n");
		exit(1);
	}

	if (!write_hyperspectral_file_to_stdout) {
		output_hyperspectral_basename = output_basename;
	}

	if ((spectral_coordinates.size() > 0) && filename.empty()) {
		fprintf(stderr, "Cannot use image coordinates as masking spectra when hyperspectral image is arriving from stdin.\n");
		exit(1);
	}

	//read image
	struct hyspex_header header;
	check_reading_errors(filename, hyperspectral_read_header(filename.c_str(), &header));
	hyperspectral_reader_t *reader = hyperspectral_create_reader(header, filename);
	float *line_data = new float[header.samples*header.bands]();

	masking_t mask_param;
	masking_err_t errcode;

	if (spectral_coordinates.size() > 0) {
		//extract masking spectra from image
		std::vector<float*> spectra;
		for (int i=0; i < spectral_coordinates.size(); i++) {
			check_reading_errors(filename, hyperspectral_line_data_float(reader, spectral_coordinates[i].first, line_data));
			spectra.push_back(new float[header.bands]());
			for (int j=0; j < header.bands; j++) {
				spectra[i][j] = line_data[j*header.samples + spectral_coordinates[i].second];
			}
		}

		errcode = masking_init(header.wlens.size(), spectra, &mask_param, sam_thresh);
		for (int i=0; i < spectra.size(); i++) {
			delete [] spectra[i];
		}
	} else if (spectral_files.size() > 0) {
		//get masking spectra from files
		errcode = masking_init(header.wlens, spectral_files, &mask_param, sam_thresh);
	} else if (!spectral_directory.empty()) {
		//get masking spectra from directory
		errcode = masking_init(header.wlens, spectral_directory, &mask_param, sam_thresh);
	} else {
		fprintf(stderr, "No masking spectra defined.\n");
		exit(1);
	}

	if (errcode != MASKING_NO_ERR) {
		fprintf(stderr, "Error in initializing masking parameters: %s\n", masking_error_message(errcode));
		exit(1);
	}

	mask_thresh_t thresh_val = masking_allocate_thresh(&mask_param, header.samples);

	hyperspectral_file_t file;

	if (write_hyperspectral_output) {
		file = hyperspectral_open(output_hyperspectral_basename, header.lines, header.bands, header.samples, header.wlens);
	}

	std::ofstream *mask_file;
	if (write_mask_to_file) {
		mask_file = new std::ofstream(output_basename + "_mask.dat");
	}

	//read image and mask
	for (int i=0; i < header.lines; i++){
		if (verbosity > 0) print_progress_bar("Mask image", i-1, i, header.lines-1);

		check_reading_errors(filename, hyperspectral_line_data_float(reader, i, line_data));

		//mask line
		masking_thresh(&mask_param, header.samples, line_data, &thresh_val);

		for (int i=0; i < header.samples; i++){
			if (write_hyperspectral_output) {
				for (int j=0; j < header.bands; j++) {
					if (!masking_pixel_belongs(&mask_param, thresh_val, i)) {
						line_data[j*header.samples + i] = 0;
					}
				}
			}

			if (write_mask_to_file) {
				*mask_file << masking_pixel_belongs(&mask_param, thresh_val, i) << " ";
			}
		}

		if (write_mask_to_file) {
			*mask_file << std::endl;
		}

		//output line
		if (write_hyperspectral_output) {
			hyperspectral_write_to_file(&file, header.samples*header.bands*sizeof(float), (char*)line_data);
		}
	}
	if (write_hyperspectral_output) {
		hyperspectral_close_write_file(&file);
	}

	if (write_mask_to_file) {
		mask_file->close();
	}
	masking_free_thresh(&thresh_val, header.samples);
	masking_free(&mask_param);
	delete [] line_data;
}
