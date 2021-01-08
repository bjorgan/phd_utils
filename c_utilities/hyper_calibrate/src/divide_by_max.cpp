#include "divide_by_max.h"
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>

void calibrate_line_by_maximum_value(struct hyspex_header header, float *line)
{
	for (int i=0; i < header.samples; i++) {
		float maxvalue = 0;
		for (int j=0; j < header.bands; j++) {
			float value = line[j*header.samples + i];
			if (value > maxvalue) {
				maxvalue = value;
			}
		}
		
		for (int j=0; j < header.bands; j++) {
			size_t index = j*header.samples + i;
			line[index] = line[index]/maxvalue;
		}
	}
}

void calibrate_by_max_spectral_value(std::string input_filename, std::string output_filename)
{
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header(input_filename.c_str(), &header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
		exit(1);
	}

	hyperspectral_reader_t *reader = hyperspectral_create_reader(header, input_filename);
	hyperspectral_file_t writer = hyperspectral_open(output_filename, header.lines, header.bands, header.samples, header.wlens);

	float *line = new float[header.samples*header.bands]();
	size_t line_size = sizeof(float)*header.samples*header.bands;

	for (int i=0; i < header.lines; i++) {
		errcode = hyperspectral_line_data_float(reader, i, line);
		if (errcode != HYPERSPECTRAL_NO_ERR) {
			fprintf(stderr, "Error reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
			exit(1);
		}

		calibrate_line_by_maximum_value(header, line);
		size_t written_size = hyperspectral_write_to_file(&writer, line_size, (char*)line);
		if (written_size < line_size) {
			fprintf(stderr, "Stream closed.\n");
			exit(1);
		}
	}

	hyperspectral_close_write_file(&writer);
}
