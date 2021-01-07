#include "generate_image.hpp"
#include <iostream>
#include "readimage.h"
#include "readimage_stdin.h"

#define CHECK_EQUAL(expected, actual) \
{ if ((expected) != (actual)) { \
	std::cerr << "Failed on line " << __FILE__ << ":" << __LINE__ << ", expected " << expected << ", got " << actual << std::endl; \
	exit(1);\
  } }

#define CHECK_INEQUAL(expected, actual) \
{ if ((expected) == (actual)) { \
	std::cerr << "Failed on line " << __FILE__ << ":" << __LINE__ << ", expected " << expected << " to be different from " << actual << std::endl; \
	exit(1);\
  } }

int main()
{
	//generate original image
	int num_lines = 0, num_bands = 0, num_samples = 0;
	std::vector<float> wlens;
	float *orig_image = generate_image(&num_lines, &num_bands, &num_samples, &wlens);

	//read header from stdin
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header_from_stdin(&header);
	CHECK_EQUAL(HYPERSPECTRAL_NO_ERR, errcode);
	CHECK_EQUAL(num_lines, header.lines);
	CHECK_EQUAL(num_samples, header.samples);
	CHECK_EQUAL(num_bands, header.bands);

	//should be equal on first read
	hyperspectral_stdin_t *data = hyperspectral_stdin_create(header, STDIN_ACCUMULATE_DATA);
	for (int i=0; i < header.lines; i++) {
		const float *input_line = (float*)hyperspectral_stdin_line_data(data, i);
		for (int j=0; j < header.samples*header.bands; j++) {
			CHECK_EQUAL(input_line[j], orig_image[i*header.samples*header.bands + j]);
		}
	}

	//should be equal on second read
	for (int i=0; i < header.lines; i++) {
		CHECK_INEQUAL(0, hyperspectral_stdin_line_data(data, i));

		const float *input_line = (float*)hyperspectral_stdin_line_data(data, i);
		for (int j=0; j < header.samples*header.bands; j++) {
			CHECK_EQUAL(input_line[j], orig_image[i*header.samples*header.bands + j]);
		}
	}

	delete [] orig_image;

	return 0;

}
