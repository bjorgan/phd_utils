#include "calibration_io.h"
#include <fstream>
#include <sstream>
#include <iostream>

const std::string NUM_BANDS_PROPERTY = "num_bands";
const std::string NUM_SAMPLES_PROPERTY = "num_samples";
const std::string ADJUSTED_SPECS_PROPERTY = "adjusted_against_specifications";
const std::string DATA_PROPERTY = "data";

const std::string START_LINE = "start_line";
const std::string END_LINE = "end_line";
const std::string START_SAMPLE = "start_sample";
const std::string END_SAMPLE = "end_sample";

void hyperspectral_calibration_info_to_file(const hyperspectral_calibration_info_t *calibration, std::string filename)
{
	std::ofstream file;
	std::ostream &output = filename.empty() ? std::cout : file;

	if (!filename.empty()) {
		file.open(filename.c_str());
	}

	output.precision(std::numeric_limits<double>::digits10 + 1);

	//header info
	output << NUM_BANDS_PROPERTY << " = " << calibration->num_bands << std::endl;
	output << NUM_SAMPLES_PROPERTY << " = " << calibration->num_samples << std::endl;
	output << ADJUSTED_SPECS_PROPERTY << " = " << calibration->adjusted << std::endl;

	//standard position info
	output << START_LINE << " = " << calibration->refl_standard_positions.start_line << std::endl;
	output << END_LINE << " = " << calibration->refl_standard_positions.end_line << std::endl;
	output << START_SAMPLE << " = " << calibration->refl_standard_positions.start_sample << std::endl;
	output << END_SAMPLE << " = " << calibration->refl_standard_positions.end_sample << std::endl;

	//calibration data
	output << DATA_PROPERTY << " = {";
	for (int i=0; i < calibration->num_bands*calibration->num_samples; i++) {
		output << calibration->calibration_array[i] << " ";
	}
	output << std::endl << "}" << std::endl;

	if (!filename.empty()) {
		file.close();
	}
}

std::string header_extract_property(std::string, std::string); //from readimage.cpp

hyperspectral_calibration_info_t hyperspectral_calibration_info_from_file(std::string filename)
{
	hyperspectral_calibration_info_t ret_info = {0};

	//read file
	std::ifstream file(filename.c_str());
	std::stringstream buffer;
	buffer << file.rdbuf();
	std::string calibration_text = buffer.str();
	file.close();
	if (!calibration_text.empty()) {
		//get properties
		std::string num_bands = header_extract_property(calibration_text, NUM_BANDS_PROPERTY);
		std::string num_samples = header_extract_property(calibration_text, NUM_SAMPLES_PROPERTY);
		std::string adjusted = header_extract_property(calibration_text, ADJUSTED_SPECS_PROPERTY);
		std::string data = header_extract_property(calibration_text, DATA_PROPERTY);

		ret_info.refl_standard_positions.start_line = atoi(header_extract_property(calibration_text, START_LINE).c_str());
		ret_info.refl_standard_positions.end_line = atoi(header_extract_property(calibration_text, END_LINE).c_str());
		ret_info.refl_standard_positions.start_sample = atoi(header_extract_property(calibration_text, START_SAMPLE).c_str());
		ret_info.refl_standard_positions.end_sample = atoi(header_extract_property(calibration_text, END_SAMPLE).c_str());

		ret_info.num_bands = atoi(num_bands.c_str());
		ret_info.num_samples = atoi(num_samples.c_str());
		ret_info.calibration_array = new double[ret_info.num_bands*ret_info.num_samples]();
		ret_info.adjusted = atoi(adjusted.c_str());

		//get calibration array
		int data_pos = data.find("{") + 1;
		std::stringstream data_stream(data.substr(data_pos));
		double value = 0;
		int i=0;
		while (data_stream >> value) {
			ret_info.calibration_array[i] = value;
			i++;
		}
	}

	return ret_info;
}

bool hyperspectral_calibration_info_valid(struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info)
{
	return (calibration_info->calibration_array != NULL) && (calibration_info->num_bands == header.bands) && ((calibration_info->num_samples == header.samples) || (calibration_info->num_samples == 1));
}
