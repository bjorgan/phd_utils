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
	std::cerr << "Failed on line " << __FILE__ << ":" << __LINE__ << ", expected " << expected << ", got " << actual << std::endl; \
	exit(1);\
  } }

#define CHECK(expression) \
{ if (!(expression)) { \
	std::cerr << "Failed on line " << __FILE__ << ":" << __LINE__ << std::endl; \
	exit(1);\
  } }

int main(int argc, char *argv[])
{
	//generate original image
	int num_lines = 0, num_bands = 0, num_samples = 0;
	std::vector<float> wlens;
	float *orig_image = generate_image(&num_lines, &num_bands, &num_samples, &wlens);

	int start_line = 0;
	int end_line = num_lines;

	if (argc >= 2) {
		start_line = atoi(argv[1]);
	}
	if (argc >= 3) {
		end_line = atoi(argv[2]);
	}

	CHECK(start_line >= 0);
	CHECK(start_line < num_lines);
	CHECK(end_line > 0);
	CHECK(end_line <= num_lines);

	//read header from stdin
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header_from_stdin(&header);
	CHECK_EQUAL(HYPERSPECTRAL_NO_ERR, errcode);
	CHECK_EQUAL(num_lines, header.lines);
	CHECK_EQUAL(num_samples, header.samples);
	CHECK_EQUAL(num_bands, header.bands);
	CHECK_EQUAL(HYPERSPECTRAL_DATATYPE_FLOAT, header.datatype);

	//should be equal on first read
	hyperspectral_stdin_t *data = hyperspectral_stdin_create(header, STDIN_NO_ACCUMULATION);
	for (int i=start_line; i < end_line; i++) {
		const float *input_line = (float*)hyperspectral_stdin_line_data(data, i);
		CHECK(input_line != NULL);
		for (int j=0; j < header.samples*header.bands; j++) {
			CHECK_EQUAL(input_line[j], orig_image[i*header.samples*header.bands + j]);
		}
	}

	//next read should yield NULL pointers on everything
	for (int i=0; i < header.lines; i++) {
		CHECK_EQUAL(0, hyperspectral_stdin_line_data(data, i));
	}

	delete [] orig_image;

	return 0;

}
