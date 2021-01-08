#include "calibration_io.h"
#include <CppUTest/TestHarness.h>
#include <iostream>
#include <unistd.h>
#include <string.h>

TEST_GROUP(CalibrationIO)
{
	char *tmp_dir;
	std::string tmp_output_base;
	std::string tmp_output_file;

	void setup() {
		char dir_template[] = "/tmp/caltest.XXXXXX";
		tmp_dir = strdup(mkdtemp(dir_template));
		tmp_output_base = std::string(tmp_dir) + "/test_cal";
		tmp_output_file = tmp_output_base;
	}

	void teardown() {
		CHECK_EQUAL(0, unlink(tmp_output_file.c_str()));
		CHECK_EQUAL(0, rmdir(tmp_dir));
		free(tmp_dir);
	}
};

TEST(CalibrationIO, testReadWrite)
{
	//test writing
	hyperspectral_calibration_info_t cal;
	cal.num_samples = 10;
	cal.num_bands = 35;
	cal.adjusted = true;
	cal.calibration_array = new double[cal.num_samples*cal.num_bands]();
	for (int i=0; i < cal.num_bands*cal.num_samples; i++) {
		cal.calibration_array[i] = i*5.6739;
	}

	hyperspectral_calibration_info_to_file(&cal, tmp_output_base);
	CHECK_EQUAL(0, access(tmp_output_file.c_str(), F_OK));

	//test reading
	hyperspectral_calibration_info_t read_cal = hyperspectral_calibration_info_from_file(tmp_output_file);

	CHECK_EQUAL(cal.num_samples, read_cal.num_samples);
	CHECK_EQUAL(cal.num_bands, read_cal.num_bands);
	CHECK_EQUAL(cal.adjusted, read_cal.adjusted);

	const double TOLERANCE = 0.00001;
	for (int i=0; i < cal.num_bands*cal.num_samples; i++) {
		DOUBLES_EQUAL(cal.calibration_array[i], read_cal.calibration_array[i], TOLERANCE);
	}

	hyperspectral_calibration_info_free(&cal);
	hyperspectral_calibration_info_free(&read_cal);
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
