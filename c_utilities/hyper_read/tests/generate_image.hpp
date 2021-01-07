#include <vector>
#include <stdio.h>

float *generate_image(int *num_lines, int *num_bands, int *num_samples, std::vector<float> *wlens)
{
	*num_lines = 20;
	*num_samples = 30;
	*num_bands = 37;
	for (int i=0; i < *num_bands; i++) {
		wlens->push_back((i + 45)*3.75);
	}

	size_t size = (*num_lines)*(*num_samples)*(*num_bands);
	float *image = new float[size]();
	for (int i=0; i < size; i++) {
		image[i] = (i + 23)*5.6;
	}

	return image;
}
