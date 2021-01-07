#include "generate_image.hpp"
#include <iostream>
#include "readimage.h"

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

int main()
{
	//generate original image
	int num_lines = 0, num_bands = 0, num_samples = 0;
	std::vector<float> wlens;
	float *orig_image = generate_image(&num_lines, &num_bands, &num_samples, &wlens);

	//read header from stdio
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header_from_stdin(&header);
	CHECK_EQUAL(HYPERSPECTRAL_NO_ERR, errcode);
	CHECK_EQUAL(num_lines, header.lines);
	CHECK_EQUAL(num_samples, header.samples);
	CHECK_EQUAL(num_bands, header.bands);

	//read image from stdio
	float *read_image = hyperspectral_alloc_float(header);
	int line_number = 0;
	while (true) {
		float *line = read_image + line_number*header.samples*header.bands;
		int num_read_elements = fread(line, sizeof(float), header.samples*header.bands, stdin);
		if (num_read_elements == 0) {
			break;
		}
		line_number++;

		if (line_number > header.lines) {
			fprintf(stderr, "Unexpected number of lines.\n");
			exit(1);
		}

		if (num_read_elements != header.samples*header.bands) {
			fprintf(stderr, "Unexpected number of read bytes.\n");
			exit(1);
		}
	}

	//check against read header and image
	for (int i=0; i < num_bands; i++) {
		CHECK_EQUAL(wlens[i], header.wlens[i]);
	}

	for (int i=0; i < num_bands*num_samples*num_lines; i++) {
		CHECK_EQUAL(orig_image[i], read_image[i]);
	}

	delete [] orig_image;
	delete [] read_image;

	return 0;

}
