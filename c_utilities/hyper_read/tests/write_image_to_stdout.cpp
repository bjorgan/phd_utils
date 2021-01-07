#include "generate_image.hpp"
#include "readimage.h"

int main() {
	int num_lines, num_samples, num_bands;
	std::vector<float> wlens;
	float *image = generate_image(&num_lines, &num_bands, &num_samples, &wlens);

	hyperspectral_file_t file = hyperspectral_open_stdout_file(num_lines, num_bands, num_samples, wlens);
	hyperspectral_write_to_file(&file, sizeof(float)*num_lines*num_bands*num_samples, (char*)image);
	hyperspectral_close_write_file(&file);

	return 0;
}
