#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include "getopt_long_helpers.h"
#include "calibration_io.h"
#include <unistd.h>

hyperspectral_calibration_info_t read_spectral_file(std::string filename, std::vector<float> *wlens);

int main(int argc, char *argv[])
{
	//filename for reflectance standard specification file (assumed to be explicitly set if non-empty)
	std::string specs_filename;

	//base output filename
	std::string out_filename;

	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"specs-filename",		required_argument,	0,	's'},
		{"output-base",			required_argument,	0,	'o'},
		{0, 0, 0, 0}
	};
	char short_options[] = "h:s:o";
	const char *option_descriptions[] = {
		"Show help",
		"Filename of specification file for what the reflectance from the reflectance standard is. It is assumed that the specification file consists of two columns, wavelengths and values between 0 and 1",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] spectral_file\n";

	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;

			case 's': //reflectance standard specification
				specs_filename = std::string(optarg);
			break;

			case 'o':
				out_filename = std::string(optarg);
			break;
		}
		if (flag == -1){
			break;
		}
	}

	//input filename
	std::string filename;
	if (optind < argc) {
		filename = std::string(argv[optind]);
	}

	//read spectra from file
	int num_samples, num_bands;
	std::vector<float> wlens;
	hyperspectral_calibration_info_t calibration_info = read_spectral_file(filename, &wlens);

	struct hyspex_header dummy_header;
	dummy_header.samples = num_samples;
	dummy_header.bands = num_bands;
	dummy_header.wlens = wlens;

	//adjust against specifications
	if (!specs_filename.empty()){
		calibration_adjust(specs_filename, dummy_header, &calibration_info);
	}

	//write calibration info to file
	hyperspectral_calibration_info_to_file(&calibration_info, out_filename);
}

hyperspectral_calibration_info_t read_spectral_file(std::string filename, std::vector<float> *wlens)
{
	std::vector<std::vector<double> > data;

	//read file/stdin input to vector
	std::ifstream file_input;
	std::istream &input = filename.empty() ? std::cin : file_input;

	if (!filename.empty()) {
		file_input.open(filename);
	}

	std::string line;
	while(!std::getline(input, line).eof()) {
		std::istringstream reader(line);
		std::vector<double> line_data;
		while (!reader.eof()) {
			double val = 0;
			reader >> val;
			if (reader.fail()) {
				break;
			}

			line_data.push_back(val);
		}
		data.push_back(line_data);
	}

	if (!filename.empty()) {
		file_input.close();
	}

	//parse data input: assuming cols [wlens data data data], with new wavelength on every row.
	int num_samples = data[0].size() - 1;
	int num_bands = data.size();
	double *ret_calibration_array = new double[num_samples*num_bands]();
	for (int i=0; i < data.size(); i++) {
		wlens->push_back(data[i][0]);
		for (int j=1; j < data[i].size(); j++) {
			ret_calibration_array[i*num_samples + j - 1] = data[i][j];
		}
	}

	//generate calibration info
	struct image_subset dummy_position = {0};
	hyperspectral_calibration_info_t calibration_info = calibration_create_info(num_samples, num_bands, ret_calibration_array, dummy_position);
	delete [] ret_calibration_array;

	return calibration_info;
}

