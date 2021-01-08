#include <hyperspectral/readimage.h>
#include <CppUTest/TestHarness.h>
#include <iostream>
#include <random>
#include <string.h>

#include "mnf.h"

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
float *global_noise_image;

/**
 * Generate test hyperspectral image. Uses LINES, SAMPLES and BANDS to define the sizes, and MEAN and STDEV to generate
 * gaussian noise. Divides the images in four patches, with the following spectra in each patch: cosine, sine, straight
 * line with positive derivative, straight line with negative derivative.
 **/
void generate_image(struct hyspex_header *output_header, float **noise_free, float **noisy) {
	output_header->lines = LINES;
	output_header->samples = SAMPLES;
	output_header->bands = BANDS;

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
	float *noise_image = hyperspectral_alloc_float(*output_header);
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(MEAN, STDEV);
	for (int i=0; i < output_header->lines*output_header->samples*output_header->bands; i++) {
		double noise = distribution(generator);
		noise_image[i] = image[i] + noise;
	}

	*noisy = noise_image;
}

///Threshold for when we consider to values to be equal
double EQUALITY_THRESHOLD = 0.0001;

/**
 * Check that MNF yields the identity transform when the number of bands to use in inverse
 * equals all bands.
 **/
TEST_GROUP(mnf_equality)
{
	struct hyspex_header header;
	float *noise_image;
	mnf_statistics_t *statistics;
	mnf_transform_t *transform;
	size_t all_samples;

	void setup() {
		//copy noisy image from global array
		header = global_header;
		noise_image = hyperspectral_alloc_float(header);
		memcpy(noise_image, global_noise_image, sizeof(float)*header.samples*header.lines*header.bands);
		all_samples = header.samples*header.lines*header.bands;

		//create and calculate statistics
		statistics = mnf_statistics_create(header);
		mnf_update_statistics(statistics, header, noise_image);

		//create and calculate mnf transform
		transform = mnf_transform_create(header, header.bands);
		CHECK_EQUAL(MNF_NO_ERR, mnf_transform_calculate(header, statistics, transform));
	}

	//check image equality (transform should be more or less the identity transform)
	void image_compare(double tolerance) {
		for (int i=0; i < all_samples; i++) {
			DOUBLES_EQUAL(global_noise_image[i], noise_image[i], tolerance);
		}
	}

	//check that the eigenvalues make sense
	void eigenvalue_compare() {
		bool all_are_equal = true;
		for (int i=1; i < header.bands; i++) {
			CHECK(transform->eigenvalues[i] >= transform->eigenvalues[i-1]);

			if (transform->eigenvalues[i] == transform->eigenvalues[i-1]) {
				all_are_equal = false;
			}
		}
		CHECK(all_are_equal);
	}

	void teardown() {
		image_compare(EQUALITY_THRESHOLD);
		mnf_free_statistics(&statistics);
		mnf_free_transform(&transform);
		delete [] noise_image;
	}
};

/**
 * Check that combining forward with inverse yields output equal to input on identity transform.
 **/
TEST(mnf_equality, combined_transforms)
{
	mnf_run_full(transform, header, noise_image);
}

/**
 * Check that forward followed by inverse yields output equal to input when transform should be the identity transform.
 **/
TEST(mnf_equality, separate_transforms)
{
	mnf_run_forward(transform, header, noise_image);
	mnf_run_inverse(transform, header, noise_image);
}

/**
 * Check that line-by-line variant yields output equal to input when transform should be the identity transform.
 **/
TEST(mnf_equality, line_by_line)
{
	//create new transform in order to let line-by-line-version update the statistics and recalculate the transform
	mnf_transform_t *transform = mnf_transform_create(header, header.bands);
	mnf_statistics_t *statistics = mnf_statistics_create(header);

	for (int i=0; i < header.lines; i++) {
		CHECK_EQUAL(MNF_NO_ERR, mnf_linebyline(transform, statistics, header, noise_image + i*header.bands*header.samples));
	}

	mnf_free_transform(&transform);
	mnf_free_statistics(&statistics);
}

const int NUM_IN_INVERSE_FOR_DENOISING = 2;

/**
 * Test general denoising capabilities. Only checking general sanity and that mnf-lbl and mnf are equal in the end.
 **/
TEST_GROUP(mnf_denoising)
{
	struct hyspex_header header;
	float *noise_image;
	mnf_statistics_t *statistics;
	mnf_transform_t *transform;
	size_t all_samples;

	void setup() {
		//copy noisy image from global array
		header = global_header;
		noise_image = hyperspectral_alloc_float(header);
		all_samples = header.samples*header.lines*header.bands;
		memcpy(noise_image, global_noise_image, sizeof(float)*all_samples);

		//create and calculate statistics
		statistics = mnf_statistics_create(header);
		mnf_update_statistics(statistics, header, noise_image);

		//create and calculate mnf transform
		transform = mnf_transform_create(header, NUM_IN_INVERSE_FOR_DENOISING);
		CHECK_EQUAL(MNF_NO_ERR, mnf_transform_calculate(header, statistics, transform));
	}

	void image_compare(double tolerance) {
		for (int i=0; i < all_samples; i++) {
			DOUBLES_EQUAL(global_image[i], noise_image[i], tolerance);
		}
	}

	void teardown() {
		image_compare(0.5); //test that images are equal-ish, difficult to verify that they give a denoised impression.
		mnf_free_statistics(&statistics);
		mnf_free_transform(&transform);
		delete [] noise_image;
	}
};

/**
 * Test denoising in general for cascading transforms.
 **/
TEST(mnf_denoising, combined_transforms)
{
	mnf_run_full(transform, header, noise_image);
}

/**
 * Test denoising in general for forward run before inverse.
 **/
TEST(mnf_denoising, separate_transforms)
{
	mnf_run_forward(transform, header, noise_image);
	mnf_run_inverse(transform, header, noise_image);
}

/**
 * Test that the last line in mnf-lbl is equal to the last line of mnf.
 **/
TEST(mnf_denoising, mnf_line_by_line_end_result_combined)
{
	float *noise_image_copy = new float[all_samples]();
	memcpy(noise_image_copy, noise_image, all_samples*sizeof(float));

	//calculate conventional transform
	mnf_run_full(transform, header, noise_image);

	//run lbl
	mnf_transform_t *transform_lbl = mnf_transform_create(header, NUM_IN_INVERSE_FOR_DENOISING);
	mnf_statistics_t *statistics_lbl = mnf_statistics_create(header);
	for (int i=0; i < header.lines; i++) {
		CHECK_EQUAL(MNF_NO_ERR, mnf_linebyline(transform_lbl, statistics_lbl, header, noise_image_copy + i*header.bands*header.samples));
	}
	mnf_free_transform(&transform_lbl);
	mnf_free_statistics(&statistics_lbl);

	//compare last line of results
	float *last_full_line = noise_image + (header.lines-1)*header.samples*header.bands;
	float *last_lbl_line = noise_image_copy + (header.lines-1)*header.samples*header.bands;
	for (int i=0; i < header.samples*header.bands; i++) {
		DOUBLES_EQUAL(last_full_line[i], last_lbl_line[i], EQUALITY_THRESHOLD);
	}
}

#include <unistd.h>
#include <string.h>

/**
 * Test MNF statistics IO.
 **/
TEST_GROUP(mnf_io)
{
	struct hyspex_header header;
	float *noise_image;
	mnf_statistics_t *statistics;
	char *tmp_dir;
	std::string basename;
	std::string bandmeans;
	std::string noisecov;
	std::string imagecov;

	void setup() {
		//create a temporary directory to which we can write statistics
		char dir_template[] = "/tmp/hypermnftest.XXXXXX";
		tmp_dir = strdup(mkdtemp(dir_template));
		CHECK(tmp_dir != NULL);
		basename = std::string(tmp_dir) + "/mnf_statistics_test";
		bandmeans = basename + MNF_BANDMEANS_FILE_POSTFIX;
		noisecov = basename + MNF_NOISECOV_FILE_POSTFIX;
		imagecov = basename + MNF_IMAGECOV_FILE_POSTFIX;

		//copy noisy image from global array
		header = global_header;
		noise_image = hyperspectral_alloc_float(header);
		size_t all_samples = header.samples*header.lines*header.bands;
		memcpy(noise_image, global_noise_image, sizeof(float)*all_samples);

		//create and calculate statistics
		statistics = mnf_statistics_create(header);
		mnf_update_statistics(statistics, header, noise_image);
	}

	bool file_exists() {
		return (access(bandmeans.c_str(), F_OK) == 0)
			&& (access(noisecov.c_str(), F_OK) == 0)
			&& (access(imagecov.c_str(), F_OK) == 0);
	}

	void teardown() {
		mnf_free_statistics(&statistics);
		delete [] noise_image;

		//remove temporary files
		CHECK_EQUAL(0, unlink(bandmeans.c_str()));
		CHECK_EQUAL(0, unlink(noisecov.c_str()));
		CHECK_EQUAL(0, unlink(imagecov.c_str()));
		CHECK_EQUAL(0, rmdir(tmp_dir));
		free(tmp_dir);
	}
};

TEST(mnf_io, write_to_file)
{
	//write to file
	CHECK(!mnf_statistics_exist(basename));
	CHECK(!file_exists());
	mnf_statistics_to_file(basename, header, statistics);
	CHECK(file_exists());
	CHECK(mnf_statistics_exist(basename));
}

#include "mnf_private.h"

void covariances_equal(struct hyspex_header header, imagestatistics_t *written_stats, imagestatistics_t *read_stats, std::string message)
{
	//first sanity check
	bool all_cov_zero = true;
	bool all_means_zero = true;
	for (int i=0; i < header.bands*header.bands; i++) {
		if (read_stats->C[i] != 0.0f) all_cov_zero = false;
	}
	CHECK_TEXT(read_stats->n != 0, std::string("Number of samples was zero, " + message).c_str());
	CHECK_TEXT(!all_cov_zero, std::string("All covariance values were zero, " + message).c_str());

	//calculate covariances
	float *cov_written = new float[header.bands*header.bands]();
	imagestatistics_get_cov(written_stats, header.bands, cov_written);
	float *cov_read = new float[header.bands*header.bands]();
	imagestatistics_get_cov(read_stats, header.bands, cov_read);

	//compare struct cov and cov gotten from get_cov, as in case of read stats this should be the same
	for (int i=0; i < header.bands*header.bands; i++) {
		DOUBLES_EQUAL_TEXT(read_stats->C[i], cov_read[i], EQUALITY_THRESHOLD, std::string("Covariance sanity, " + message).c_str());
	}

	//compare covariances
	for (int i=0; i < header.bands*header.bands; i++) {
		DOUBLES_EQUAL_TEXT(cov_written[i], cov_read[i], EQUALITY_THRESHOLD, std::string("Covariance equality between written and read statistics, " + message).c_str());
	}

	//compare ones_samples
	for (int i=0; i < header.samples; i++) {
		DOUBLES_EQUAL_TEXT(written_stats->ones_samples[i], read_stats->ones_samples[i], EQUALITY_THRESHOLD, message.c_str());
	}
}

void means_equal(struct hyspex_header header, imagestatistics_t *written_stats, imagestatistics_t *read_stats, std::string message)
{
	//calculate means
	float *means_written = new float[header.bands]();
	imagestatistics_get_means(written_stats, header.bands, means_written);
	float *means_read = new float[header.bands]();
	imagestatistics_get_means(read_stats, header.bands, means_read);

	//compare means
	for (int i=0; i < header.bands; i++) {
		DOUBLES_EQUAL_TEXT(means_written[i], means_read[i], EQUALITY_THRESHOLD, std::string("Mean equality between written and read statistics, " + message).c_str());
	}
}

void transforms_equal(struct hyspex_header header, mnf_statistics_t *written_statistics, mnf_statistics_t *read_statistics)
{
	mnf_transform_t *written_transform = mnf_transform_create(header, header.bands/3);
	CHECK_EQUAL(MNF_NO_ERR, mnf_transform_calculate(header, written_statistics, written_transform));

	mnf_transform_t *read_transform = mnf_transform_create(header, header.bands/3);
	CHECK_EQUAL(MNF_NO_ERR, mnf_transform_calculate(header, read_statistics, read_transform));

	for (int i=0; i < header.bands*header.bands; i++) {
		DOUBLES_EQUAL_TEXT(written_transform->forward_transform[i], read_transform->forward_transform[i], EQUALITY_THRESHOLD, "Forward transforms");
		DOUBLES_EQUAL_TEXT(written_transform->inverse_transform[i], read_transform->inverse_transform[i], EQUALITY_THRESHOLD, "Inverse transforms");
	}

	mnf_free_transform(&written_transform);
	mnf_free_transform(&read_transform);
}

TEST(mnf_io, read_from_file)
{
	//write to file
	CHECK(!file_exists());
	mnf_statistics_to_file(basename, header, statistics);
	CHECK(file_exists());

	//read from file
	mnf_statistics_t *read_statistics = mnf_statistics_create(header);
	CHECK_EQUAL(MNF_NO_ERR, mnf_statistics_from_file(basename, header, read_statistics));

	//check that statistics are equal
	covariances_equal(header, image_statistics(statistics), image_statistics(read_statistics), "image statistics");
	means_equal(header, image_statistics(statistics), image_statistics(read_statistics), "image statistics");

	struct hyspex_header noise_header = header;
	noise_header.samples = noise_header.samples-1;
	covariances_equal(noise_header, noise_statistics(statistics), noise_statistics(read_statistics), "noise statistics");

	//test that resulting transforms are equal
	transforms_equal(header, statistics, read_statistics);
}

#include <vector>
#include <CppUTest/CommandLineTestRunner.h>

int main(int argc, char** argv)
{
	generate_image(&global_header, &global_image, &global_noise_image);

	MemoryLeakWarningPlugin::turnOffNewDeleteOverloads();

	std::vector<const char*> args(argv, argv + argc); // Insert all arguments
	args.push_back("-v"); // Set verbose mode
	args.push_back("-c"); // Set color output (OPTIONAL)

	return RUN_ALL_TESTS(args.size(), &args[0]);
}
