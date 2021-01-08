#include <hyperspectral/readimage.h>
#include <CppUTest/TestHarness.h>
#include <iostream>
#include <random>
#include <string.h>
#include "compression.h"
#include <hyperspectral/mnf.h>

/**
 * Size of hyperspectral test image.
 **/
const int LINES = 139;
const int SAMPLES = 145;
const int BANDS = 21;

/**
 * Default noise statistics.
 **/
const double MEAN = 0;
const double STDEV = sqrt(0.01);

/**
 * Original images, are copied into local images in each test.
 **/
struct hyspex_header global_header;
float *global_image;
float *global_noisy_image;

/**
 * Generate test hyperspectral image. Uses LINES, SAMPLES and BANDS to define the sizes, and MEAN and STDEV to generate
 * gaussian noise. Divides the images in four patches, with the following spectra in each patch: cosine, sine, straight
 * line with positive derivative, straight line with negative derivative.
 **/
void generate_image(struct hyspex_header *output_header, float **noise_free, float **noisy) {
	output_header->lines = LINES;
	output_header->samples = SAMPLES;
	output_header->bands = BANDS;
	output_header->datatype = HYPERSPECTRAL_DATATYPE_FLOAT;

	for (int i=0; i < output_header->bands; i++) {
		output_header->wlens.push_back(i*3.7);
	}

	float *image = hyperspectral_alloc_float(*output_header);

	//generate image
	for (int i=0; i < output_header->lines; i++) {
		for (int k=0; k < output_header->bands; k++) {
			for (int j=0; j < output_header->samples; j++) {
				int index = i*output_header->samples*output_header->bands + k*output_header->samples + j;
				if ((i < LINES/2) && (j < SAMPLES/2)) {
					image[index] = k/(1.0f*output_header->bands);
				} else if ((i < LINES/2) && (j >= SAMPLES/2)) {
					image[index] = powf(cos(0.1*k), 2);
				} else if ((i >= LINES/2) && (j < SAMPLES/2)) {
					image[index] = powf(sin(0.1*k), 2);
				} else {
					image[index] = (output_header->bands - k)/(1.0*output_header->bands);
				}
			}
		}
	}
	*noise_free = image;

	//add noise to image
	float *noisy_image = hyperspectral_alloc_float(*output_header);
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(MEAN, STDEV);
	for (int i=0; i < output_header->lines*output_header->samples*output_header->bands; i++) {
		double noise = distribution(generator);
		noisy_image[i] = image[i] + noise;
	}

	*noisy = noisy_image;
}

///Threshold for when we consider to values to be equal
double EQUALITY_THRESHOLD = 0.0001;

const int NUM_BANDS_IN_INVERSE = 2;

TEST_GROUP(mnf_compression)
{
	struct hyspex_header header;
	float *noisy_image;
	float *mnf_denoised_image;
	size_t bytes_in_image;
	hyperspectral_compression_t *compressed_image;
	mnf_transform_t *transform;

	void setup() {
		//copy noisy image from global array
		header = global_header;
		noisy_image = hyperspectral_alloc_float(header);
		bytes_in_image = sizeof(float)*header.samples*header.lines*header.bands;
		memcpy(noisy_image, global_noisy_image, bytes_in_image);

		//create compression
		compressed_image = hyperspectral_create_compression(header, (char*)noisy_image, NUM_BANDS_IN_INVERSE);
		delete [] noisy_image;

		//create MNF denoised version of the image
		mnf_denoised_image = hyperspectral_alloc_float(header);
		memcpy(mnf_denoised_image, global_noisy_image, bytes_in_image);
		transform = mnf_transform_create(header, NUM_BANDS_IN_INVERSE);
		mnf_statistics_t *statistics = mnf_statistics_create(header);
		mnf_update_statistics(statistics, header, mnf_denoised_image);
		mnf_transform_calculate(header, statistics, transform);
		mnf_run_full(transform, header, mnf_denoised_image);
	}

	void teardown() {
		delete [] mnf_denoised_image;
		hyperspectral_free_compression(&compressed_image);
	}
};

int index(const struct hyspex_header &header, int line, int band, int sample)
{
	return line*header.samples*header.bands + band*header.samples + sample;
}

TEST(mnf_compression, spectrum)
{
	float *spectrum = new float[header.bands]();
	float *line = new float[header.bands*header.samples]();

	for (int i=0; i < header.lines; i++) {
		hyperspectral_line(i, compressed_image, line);
		for (int j=0; j < header.samples; j++) {
			hyperspectral_pixel(i, j, compressed_image, spectrum);
			for (int k=0; k < header.bands; k++) {
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], line[k*header.samples + j], EQUALITY_THRESHOLD);
				DOUBLES_EQUAL(line[k*header.samples + j], spectrum[k], EQUALITY_THRESHOLD);
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], spectrum[k], EQUALITY_THRESHOLD);
			}
		}
	}
	delete [] spectrum;
}

TEST(mnf_compression, line)
{
	float *line = new float[header.samples*header.bands]();
	for (int i=0; i < header.lines; i++) {
		hyperspectral_line(i, compressed_image, line);
		for (int j=0; j < header.bands*header.samples; j++) {
			DOUBLES_EQUAL(mnf_denoised_image[i*header.samples*header.bands + j], line[j], EQUALITY_THRESHOLD);
		}
	}
	delete [] line;
}

TEST(mnf_compression, band)
{
	float *band = new float[header.samples*header.lines]();
	for (int k=0; k < header.bands; k++) {
		hyperspectral_band(k, compressed_image, band);
		for (int i=0; i < header.lines; i++) {
			for (int j=0; j < header.samples; j++) {
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], band[i*header.samples + j], EQUALITY_THRESHOLD);
			}
		}
	}
	delete [] band;
}

TEST(mnf_compression, orthonormality)
{
	float *inverse = compressed_image->transform->inverse_transform;
	int spectra = compressed_image->num_bands_in_inverse;

	//inverse transform
	for (int i=0; i < spectra; i++) {
		for (int j=0; j < spectra; j++) {
			float dot_product_forward = 0;
			float dot_product_inverse = 0;
			for (int k=0; k < header.bands; k++) {
				int ind_1 = i*header.bands + k;
				int ind_2 = j*header.bands + k;
				dot_product_inverse += inverse[ind_1]*inverse[ind_2];
			}

			if (i == j) {
				DOUBLES_EQUAL(1.0f, dot_product_inverse, EQUALITY_THRESHOLD);
			} else {
				DOUBLES_EQUAL(0.0f, dot_product_inverse, EQUALITY_THRESHOLD);
			}
		}
	}

	//inverse transform, premultiplied
	inverse = compressed_image->transform->inverse_transform_with_band_selection;
	for (int i=0; i < spectra; i++) {
		for (int j=0; j < spectra; j++) {
			float dot_product_forward = 0;
			float dot_product_inverse = 0;
			for (int k=0; k < header.bands; k++) {
				int ind_1 = k*header.bands + i;
				int ind_2 = k*header.bands + j;
				dot_product_inverse += inverse[ind_1]*inverse[ind_2];
			}

			if (i == j) {
				DOUBLES_EQUAL(1.0f, dot_product_inverse, EQUALITY_THRESHOLD);
			} else {
				DOUBLES_EQUAL(0.0f, dot_product_inverse, EQUALITY_THRESHOLD);
			}
		}
	}
}

double dot_product(std::vector<float> vec1, std::vector<float> vec2)
{
	double retval = 0.0;
	for (int i=0; i < vec1.size(); i++) {
		retval += vec1[i]*vec2[i];
	}
	return retval;
}

std::vector<float> spectrum(int line, int sample, const struct hyspex_header &header, float *image)
{
	std::vector<float> retvec(header.bands, 0);
	for (int i=0; i < header.bands; i++) {
		retvec[i] = image[line*header.samples*header.bands + i*header.samples + sample];
	}
	return retvec;
}

TEST(mnf_compression, spectrum_convenience_function)
{
	int line = 3, sample = 2;
	std::vector<float> spec = spectrum(line, sample, header, mnf_denoised_image);
	for (int i=0; i < header.bands; i++) {
		CHECK_EQUAL(spec[i], mnf_denoised_image[line*header.samples*header.bands + i*header.samples + sample]);
	}
}

TEST(mnf_compression, dot_product)
{
	CHECK_EQUAL(0, dot_product(std::vector<float>(2, 0), std::vector<float>(2, 1)));
	CHECK_EQUAL(4, dot_product(std::vector<float>(2, 2), std::vector<float>(2, 1)));

	int samples = header.lines*header.samples/16;
	CHECK(samples > 2);
	for (int i=0; i < samples; i++) {
		int line_1 = i/header.samples;
		int sample_1 = i % header.samples;
		CHECK(line_1 < header.lines);
		CHECK(sample_1 < header.samples);
		std::vector<float> spectrum_1 = spectrum(line_1, sample_1, header, mnf_denoised_image);
		for (int j=0; j < samples; j++) {
			int line_2 = j/header.samples;
			int sample_2 = j % header.samples;

			CHECK(line_2 < header.lines);
			CHECK(sample_2 < header.samples);

			std::vector<float> spectrum_2 = spectrum(line_2, sample_2, header, mnf_denoised_image);
			DOUBLES_EQUAL(dot_product(spectrum_1, spectrum_2), hyperspectral_dot_product(compressed_image, coordinate(line_1, sample_1), coordinate(line_2, sample_2)), EQUALITY_THRESHOLD);
		}
	}

}

double euclidean_distance_squared(std::vector<float> vec1, std::vector<float> vec2)
{
	double retval = 0.0;
	for (int i=0; i < vec1.size(); i++) {
		retval += pow(vec1[i] - vec2[i], 2);
	}
	return retval;
}

TEST(mnf_compression, euclidean_distance)
{
	CHECK_EQUAL(2, euclidean_distance_squared(std::vector<float>(2, 0), std::vector<float>(2, 1)));
	CHECK_EQUAL(2, euclidean_distance_squared(std::vector<float>(2, 2), std::vector<float>(2, 1)));

	int samples = header.lines*header.samples/16;
	CHECK(samples > 2);
	for (int i=0; i < samples; i++) {
		int line_1 = i/header.samples;
		int sample_1 = i % header.samples;
		CHECK(line_1 < header.lines);
		CHECK(sample_1 < header.samples);
		std::vector<float> spectrum_1 = spectrum(line_1, sample_1, header, mnf_denoised_image);
		for (int j=0; j < samples; j++) {
			int line_2 = j/header.samples;
			int sample_2 = j % header.samples;

			CHECK(line_2 < header.lines);
			CHECK(sample_2 < header.samples);

			std::vector<float> spectrum_2 = spectrum(line_2, sample_2, header, mnf_denoised_image);
			DOUBLES_EQUAL(euclidean_distance_squared(spectrum_1, spectrum_2), hyperspectral_euclidean_distance_squared(compressed_image, coordinate(line_1, sample_1), coordinate(line_2, sample_2)), EQUALITY_THRESHOLD);
		}
	}

}

#include <unistd.h>
#include <string.h>

TEST_GROUP(mnf_compression_io)
{
	struct hyspex_header header;
	float *noisy_image;
	float *mnf_denoised_image;
	size_t bytes_in_image;
	hyperspectral_compression_t *compressed_image;
	mnf_transform_t *transform;
	char *tmp_dir;
	std::string output_filename;

	void setup() {
		//create a temporary directory to which we can write statistics
		char dir_template[] = "/tmp/comprtest.XXXXXX";
		tmp_dir = strdup(mkdtemp(dir_template));
		CHECK(tmp_dir != NULL);
		output_filename = std::string(tmp_dir) + "/mnf_compression_test";

		//copy noisy image from global array
		header = global_header;
		noisy_image = hyperspectral_alloc_float(header);
		bytes_in_image = sizeof(float)*header.samples*header.lines*header.bands;
		memcpy(noisy_image, global_noisy_image, bytes_in_image);

		//create compression
		compressed_image = hyperspectral_create_compression(header, (char*)noisy_image, NUM_BANDS_IN_INVERSE);
		delete [] noisy_image;

		//create MNF denoised version of the image
		mnf_denoised_image = hyperspectral_alloc_float(header);
		memcpy(mnf_denoised_image, global_noisy_image, bytes_in_image);
		transform = mnf_transform_create(header, NUM_BANDS_IN_INVERSE);
		mnf_statistics_t *statistics = mnf_statistics_create(header);
		mnf_update_statistics(statistics, header, mnf_denoised_image);
		mnf_transform_calculate(header, statistics, transform);
		mnf_run_full(transform, header, mnf_denoised_image);
	}

	void teardown() {
		delete [] mnf_denoised_image;
		hyperspectral_free_compression(&compressed_image);
		CHECK_EQUAL(0, unlink(output_filename.c_str()));
		CHECK_EQUAL(0, rmdir(tmp_dir));
		free(tmp_dir);
	}
};

TEST(mnf_compression_io, compression_to_from_file)
{
	CHECK_TEXT(access(output_filename.c_str(), F_OK) != 0, "Output file exists too soon");
	hyperspectral_compression_to_file(output_filename, compressed_image);
	CHECK_TEXT(access(output_filename.c_str(), F_OK) == 0, "Output file was not written");
	hyperspectral_compression_t *read_compressed_image = hyperspectral_compression_from_file(output_filename);
	CHECK_TEXT(read_compressed_image != NULL, "Read image is NULL");

	float *spectrum = new float[header.bands]();
	float *line = new float[header.bands*header.samples]();

	for (int i=0; i < header.lines; i++) {
		hyperspectral_line(i, read_compressed_image, line);
		for (int j=0; j < header.samples; j++) {
			hyperspectral_pixel(i, j, read_compressed_image, spectrum);
			for (int k=0; k < header.bands; k++) {
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], line[k*header.samples + j], EQUALITY_THRESHOLD);
				DOUBLES_EQUAL(line[k*header.samples + j], spectrum[k], EQUALITY_THRESHOLD);
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], spectrum[k], EQUALITY_THRESHOLD);
			}
		}
	}
	delete [] spectrum;

	for (int i=0; i < header.lines; i++) {
		hyperspectral_line(i, read_compressed_image, line);
		for (int j=0; j < header.bands*header.samples; j++) {
			DOUBLES_EQUAL(mnf_denoised_image[i*header.samples*header.bands + j], line[j], EQUALITY_THRESHOLD);
		}
	}
	delete [] line;

	float *band = new float[header.samples*header.lines]();
	for (int k=0; k < header.bands; k++) {
		hyperspectral_band(k, read_compressed_image, band);
		for (int i=0; i < header.lines; i++) {
			for (int j=0; j < header.samples; j++) {
				DOUBLES_EQUAL(mnf_denoised_image[index(header, i, k, j)], band[i*header.samples + j], EQUALITY_THRESHOLD);
			}
		}
	}
	delete [] band;
}

void mnf_transform_to_file(std::string filename, struct hyspex_header header, mnf_transform_t *transform);
void mnf_transform_from_file(std::string filename, struct hyspex_header header, mnf_transform_t *transform);

TEST(mnf_compression_io, mnf_transform_to_from_file)
{
	//write transform to file
	CHECK_TEXT(access(output_filename.c_str(), F_OK) != 0, "Output file exists too soon");
	mnf_transform_to_file(output_filename, header, transform);
	CHECK_TEXT(access(output_filename.c_str(), F_OK) == 0, "Output file was not written");

	//read transform from file
	mnf_transform_t *read_transform = mnf_transform_create(header, NUM_BANDS_IN_INVERSE);
	mnf_transform_from_file(output_filename, header, read_transform);

	for (int i=0; i < header.bands*header.bands; i++) {
		DOUBLES_EQUAL_TEXT(transform->forward_transform[i], read_transform->forward_transform[i], EQUALITY_THRESHOLD, "forward");
		DOUBLES_EQUAL_TEXT(transform->inverse_transform[i], read_transform->inverse_transform[i], EQUALITY_THRESHOLD, "inverse");
		DOUBLES_EQUAL_TEXT(transform->inverse_transform_with_band_selection[i], read_transform->inverse_transform_with_band_selection[i], EQUALITY_THRESHOLD, "band selection");
		DOUBLES_EQUAL_TEXT(transform->denoising_transform[i], read_transform->denoising_transform[i], EQUALITY_THRESHOLD, "denoising transform");
	}
	for (int i=0; i < header.bands; i++) {
		DOUBLES_EQUAL_TEXT(transform->eigenvalues[i], read_transform->eigenvalues[i], EQUALITY_THRESHOLD, "eigenvalues");
		DOUBLES_EQUAL_TEXT(transform->means[i], read_transform->means[i], EQUALITY_THRESHOLD, "means");
	}
}

TEST(mnf_compression_io, is_compressed)
{
	CHECK(access(output_filename.c_str(), F_OK) != 0);
	hyperspectral_compression_to_file(output_filename, compressed_image);
	CHECK(access(output_filename.c_str(), F_OK) == 0);

	CHECK(access((output_filename + ".hdr").c_str(), F_OK) != 0);
	hyperspectral_write_header(output_filename, header);
	CHECK(access((output_filename + ".hdr").c_str(), F_OK) == 0);

	CHECK(!hyperspectral_file_is_compressed(output_filename + ".hdr"));
	CHECK_EQUAL(unlink((output_filename + ".hdr").c_str()), 0);

	CHECK(hyperspectral_file_is_compressed(output_filename));
}



#include <vector>
#include <CppUTest/CommandLineTestRunner.h>

int main(int argc, char** argv)
{
	generate_image(&global_header, &global_image, &global_noisy_image);

	MemoryLeakWarningPlugin::turnOffNewDeleteOverloads();

	std::vector<const char*> args(argv, argv + argc); // Insert all arguments
	args.push_back("-v"); // Set verbose mode
	args.push_back("-c"); // Set color output (OPTIONAL)

	return RUN_ALL_TESTS(args.size(), &args[0]);
}
