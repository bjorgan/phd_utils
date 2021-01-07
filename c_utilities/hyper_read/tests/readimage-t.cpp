#include "readimage.cpp"
#include <CppUTest/TestHarness.h>
#include <iostream>

#define TEST_DATA_FILENAME "data/test_image.img"
#define TEST_DATA_MULTILINE_WAVELENGTH_HEADER "data/test_image_multiline_header.hdr"
#define TEST_DATA_MULTILINE_WAVELENGTH_HEADER_SPECIM "data/test_image_multiline_header_2.hdr"

TEST_GROUP(ReadImage)
{
	struct hyspex_header header;
	float *image;

	void setup() {
		image = NULL;
		CHECK_EQUAL(hyperspectral_read_header(TEST_DATA_FILENAME, &header), HYPERSPECTRAL_NO_ERR);
	}

	void teardown() {
		if (image != NULL) {
			delete [] image;
		}
	}
};

TEST(ReadImage, multilineWavelengthFieldInHeaderGetsReadCorrectly)
{
	//test with self-constructed multiline header
	struct hyspex_header multiline_header;
	CHECK_EQUAL(hyperspectral_read_header(TEST_DATA_MULTILINE_WAVELENGTH_HEADER, &multiline_header), HYPERSPECTRAL_NO_ERR);
	CHECK_EQUAL(header.samples, multiline_header.samples);
	CHECK_EQUAL(header.bands, multiline_header.bands);
	CHECK_EQUAL(header.offset, multiline_header.offset);
	CHECK_EQUAL(header.wlens.size(), multiline_header.wlens.size());
	for (int i=0; i < header.bands; i++) {
		CHECK_EQUAL(header.wlens[i], multiline_header.wlens[i]);
	}
}

void compare_image(struct hyspex_header header, const float *image, struct image_subset subset);

TEST(ReadImage, fullImage)
{
	struct image_subset full_subset = hyperspectral_generate_subset(header);
	image = hyperspectral_alloc_float(header);
	CHECK_EQUAL(hyperspectral_read_image(TEST_DATA_FILENAME, &header, image), HYPERSPECTRAL_NO_ERR);
	compare_image(header, image, full_subset);
}

TEST(ReadImage, subsetImage)
{
	struct image_subset range_subset = hyperspectral_generate_subset(header, 0, 2, 1, 3, 4, 13);
	image = hyperspectral_alloc_float(range_subset);
	CHECK_EQUAL(hyperspectral_read_image(TEST_DATA_FILENAME, &header, range_subset, image), HYPERSPECTRAL_NO_ERR);
	compare_image(header, image, range_subset);
}

TEST(ReadImage, lazyRead)
{
	hyperspectral_image_t hyperspectral_data;
	CHECK_EQUAL(hyperspectral_read_image(TEST_DATA_FILENAME, header, &hyperspectral_data), HYPERSPECTRAL_NO_ERR);
	struct image_subset subset = hyperspectral_generate_subset(header);

	const float *data = hyperspectral_float(&hyperspectral_data);
	compare_image(header, data, subset);

	hyperspectral_free_data(&hyperspectral_data);
}

TEST(ReadImage, lazyReadSubset)
{
	hyperspectral_image_t hyperspectral_data;
	CHECK_EQUAL(hyperspectral_read_image(TEST_DATA_FILENAME, header, &hyperspectral_data), HYPERSPECTRAL_NO_ERR);
	struct image_subset subset = hyperspectral_generate_subset(header);
	subset.start_line = 1;
	compare_image(header, hyperspectral_float_line(&hyperspectral_data, 1), subset);
	hyperspectral_free_data(&hyperspectral_data);
}

#define HYPER_FILENAME "test_hyperspectral_write"
#include <unistd.h>
#include <string.h>

TEST_GROUP(WriteImage)
{
	float *image;
	int num_lines, num_samples, num_bands;
	std::vector<float> wlens;
	std::string basename;
	std::string hdr_file;
	std::string img_file;
	char *tmp_dir;

	void setup() {
		//create a temporary directory to which we can write images
		char dir_template[] = "/tmp/hypertest.XXXXXX";
		tmp_dir = strdup(mkdtemp(dir_template));
		CHECK(tmp_dir != NULL);

		basename = std::string(tmp_dir) + "/" + std::string(HYPER_FILENAME);
		hdr_file = basename + ".hdr";
		img_file = basename + ".img";

		//create dummy image
		num_lines = 10;
		num_samples = 20;
		num_bands = 6;
		image = new float[num_lines*num_samples*num_bands]();
		for (int i=0; i < num_bands; i++) {
			wlens.push_back(i);
			for (int j=0; j < num_lines; j++) {
				for (int k=0; k < num_samples; k++) {
					int index = j*num_samples*num_bands + i*num_samples + k;
					image[index] = index + 4;
				}
			}
		}
	}

	void teardown() {
		//clean up temporary files and directories
		CHECK_EQUAL(unlink(img_file.c_str()), 0);
		CHECK_EQUAL(unlink(hdr_file.c_str()), 0);
		CHECK_EQUAL(rmdir(tmp_dir), 0);
		free(tmp_dir);
		delete [] image;
	}

	void check_hyperspectral_files() {
		struct hyspex_header header;
		CHECK_EQUAL(hyperspectral_read_header(img_file.c_str(), &header), HYPERSPECTRAL_NO_ERR);
		CHECK_EQUAL(header.samples, num_samples);
		CHECK_EQUAL(header.lines, num_lines);
		CHECK_EQUAL(header.bands, num_bands);
		float *read_image = new float[num_lines*num_bands*num_samples]();
		CHECK_EQUAL(hyperspectral_read_image(img_file.c_str(), &header, read_image), HYPERSPECTRAL_NO_ERR);
		for (int i=0; i < num_samples*num_lines*num_bands; i++) {
			CHECK_EQUAL(read_image[i], image[i]);
		}
	}
};

TEST(WriteImage, hyperSpectralWriteHeaderAndHyperSpectralWriteImage)
{
	CHECK(access(img_file.c_str(), F_OK ) == -1);
	CHECK(access(hdr_file.c_str(), F_OK ) == -1);

	hyperspectral_write_header(basename.c_str(), num_bands, num_samples, num_lines, wlens);
	hyperspectral_write_image(basename.c_str(), num_bands, num_samples, num_lines, image);
	CHECK(access(img_file.c_str(), F_OK ) != -1);
	CHECK(access(hdr_file.c_str(), F_OK ) != -1);

	check_hyperspectral_files();
}

TEST(WriteImage, hyperSpectralOpenWriteFileAndHyperSpectralWriteToFileAndHyperSpectralCloseWriteFile)
{
	CHECK(access(img_file.c_str(), F_OK ) == -1);
	CHECK(access(hdr_file.c_str(), F_OK ) == -1);

	hyperspectral_file_t file = hyperspectral_open_write_file(basename, num_bands, num_samples, wlens);
	CHECK(access(img_file.c_str(), F_OK ) != -1);

	hyperspectral_write_to_file(&file, sizeof(float)*num_bands*num_lines*num_samples, (char*)image);
	hyperspectral_close_write_file(&file);
	CHECK(access(hdr_file.c_str(), F_OK ) != -1);

	check_hyperspectral_files();
}

TEST_GROUP(WriteHeader)
{
};

std::string hyperspectral_header_to_string(struct hyspex_header header);

TEST(WriteHeader, hyperSpectralHeaderToString)
{
	struct hyspex_header header;
	header.interleave = BIL_INTERLEAVE;
	header.samples = 56;
	header.bands = 2;
	header.lines = 89;
	header.offset = 13;
	header.wlens.push_back(57);
	header.wlens.push_back(89);
	header.default_bands.push_back(0);
	header.default_bands.push_back(1);
	header.datatype = HYPERSPECTRAL_DATATYPE_FLOAT;

	std::string header_string = hyperspectral_header_to_string(header);
	std::cerr << header_string << std::endl;

	struct hyspex_header read_header;
	CHECK_EQUAL(HYPERSPECTRAL_NO_ERR, hyperspectral_header_from_string(header_string, &read_header));
	CHECK_EQUAL(header.interleave, read_header.interleave);
	CHECK_EQUAL(header.samples, read_header.samples);
	CHECK_EQUAL(header.bands, read_header.bands);
	CHECK_EQUAL(header.lines, read_header.lines);
	CHECK_EQUAL(header.offset, read_header.offset);
	CHECK_EQUAL(header.wlens.size(), read_header.wlens.size());
	for (int i=0; i < header.wlens.size(); i++) {
		CHECK_EQUAL(header.wlens[i], read_header.wlens[i]);
	}
	for (int i=0; i < header.default_bands.size(); i++) {
		CHECK_EQUAL(header.default_bands[i], read_header.default_bands[i]);
	}
	CHECK_EQUAL(header.datatype, read_header.datatype);
}

TEST_GROUP(ImageStdOut)
{
	float *image;
	int num_lines, num_samples, num_bands;
	std::vector<float> wlens;
	std::string stdout_filename;
	char *tmp_dir;

	void setup() {
		//create a temporary directory to which we can write images
		char dir_template[] = "/tmp/hypertest.XXXXXX";
		tmp_dir = strdup(mkdtemp(dir_template));
		CHECK(tmp_dir != NULL);
		stdout_filename = std::string(tmp_dir) + "/stdout_output.dat";

		//create dummy image
		num_lines = 10;
		num_samples = 20;
		num_bands = 6;
		image = new float[num_lines*num_samples*num_bands]();
		for (int i=0; i < num_bands; i++) {
			wlens.push_back(i);
			for (int j=0; j < num_lines; j++) {
				for (int k=0; k < num_samples; k++) {
					int index = j*num_samples*num_bands + i*num_samples + k;
					image[index] = index + 4;
				}
			}
		}
	}

	void teardown() {
		//clean up temporary files and directories
		CHECK_EQUAL(unlink(stdout_filename.c_str()), 0);
		CHECK_EQUAL(rmdir(tmp_dir), 0);
		free(tmp_dir);
		delete [] image;
	}
};

TEST(ImageStdOut, hyperspectralOpenStdoutFile)
{
	//test hyperspectral writing to stdout by writing to a specific file stream
	CHECK(access(stdout_filename.c_str(), F_OK) == -1);

	FILE *fp = fopen(stdout_filename.c_str(), "wb");
	hyperspectral_file_t file = hyperspectral_file_from_stream(fp, num_lines, num_bands, num_samples, wlens, WRITE_HEADER_TO_IMAGE_STREAM);
	
	CHECK(access(stdout_filename.c_str(), F_OK) == 0);
	
	hyperspectral_write_to_file(&file, sizeof(float)*num_bands*num_lines*num_samples, (char*)image);
	hyperspectral_close_write_file(&file);
	fclose(fp);

	//check that header was written correctly
	fp = fopen(stdout_filename.c_str(), "rb");
	struct hyspex_header header;
	CHECK_EQUAL(hyperspectral_read_header_from_binary_stream(fp, &header), HYPERSPECTRAL_NO_ERR);
	CHECK_EQUAL(num_samples, header.samples);
	CHECK_EQUAL(num_lines, header.lines);
	CHECK_EQUAL(num_bands, header.bands);
	for (int i=0; i < num_bands; i++) {
		CHECK_EQUAL(header.wlens[i], wlens[i]);
	}

	//check that file can be read correctly
	float *hyper_image = new float[num_lines*num_samples*num_bands]();
	size_t num_elements = fread((char*)hyper_image, sizeof(float), num_lines*num_samples*num_bands, fp);
	CHECK_EQUAL(num_lines*num_samples*num_bands, num_elements);
	for (int i=0; i < num_lines*num_samples*num_bands; i++) {
		CHECK_EQUAL(hyper_image[i], image[i]);
	}
}


#include <vector>
#include <CppUTest/CommandLineTestRunner.h>

int main(int argc, char** argv)
{
	MemoryLeakWarningPlugin::turnOffNewDeleteOverloads();

	std::vector<const char*> args(argv, argv + argc); // Insert all arguments
	args.push_back("-v"); // Set verbose mode
	args.push_back("-c"); // Set color output (OPTIONAL)

	return RUN_ALL_TESTS(args.size(), &args[0]);
}

void compare_image(struct hyspex_header header, const float *image, struct image_subset subset)
{
	int num_samples, num_lines, num_bands;
	hyperspectral_get_size(subset, &num_lines, &num_bands, &num_samples);

	for (int i=0; i < num_lines; i++) {
		for (int k=0; k < num_bands; k++) {
			for (int j=0; j < num_samples; j++) {
				int orig_band = k + subset.start_band;
				int orig_line = i + subset.start_line;
				int orig_sample = j + subset.start_sample;
				int orig_value = 1 + orig_band*(header.lines*header.samples) + orig_line*header.samples + orig_sample;
				int value = image[i*num_samples*num_bands + k*num_samples + j];
				CHECK_EQUAL(orig_value, value);
			}
		}
	}
}
