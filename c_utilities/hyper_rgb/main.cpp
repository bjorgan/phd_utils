#include <hyperspectral/readimage.h>
#include <iostream>
#include <getopt.h>
#include <string>
#include <cstring>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <libgen.h>
#include <unistd.h>
#include "getopt_long_helpers.h"
#include "string_helpers.h"
#include "progress_bar.h"

#define OPT_IGNORE_EXISTING 207
#define OPT_RESCALE_WITH_STATISTICS 208
#define OPT_RESCALE_TO_MAXMIN 209

std::string remove_file_ending(std::string filename); //from libhyperread
size_t hyperspectral_element_size(struct hyspex_header header); //from libhyperread

std::vector<cv::Mat> hyperspectral_to_rgb(hyperspectral_image_t *image, std::vector<int> bands, int verbosity);

int main(int argc, char *argv[])
{
	std::string output_directory;
	std::string rgb_bands_str;
	bool ignore_existing_rgb_conversions = false;
	int verbosity_level = 0;
	bool rescale_with_statistics = false;
	bool rescale_to_maxmin = false;

	//command line options
	struct option long_options[] = {
		{"output-directory",		required_argument,	0,	'o'},
		{"rgb-bands",			required_argument,	0,	'b'},
		{"help",			no_argument,		0,	'h'},
		{"verbose",			no_argument,		0,	'v'},
		{"ignore-existing",		no_argument,		0,	OPT_IGNORE_EXISTING},
		{"rescale-with-statistics",	no_argument,		0,	OPT_RESCALE_WITH_STATISTICS},
		{"rescale-to-maxmin",		no_argument,		0,	OPT_RESCALE_TO_MAXMIN},
		{0, 0, 0, 0}
	};
	char short_options[] = "o:b:h:v";

	//command line option instructions
	const char *option_descriptions[] = {
		"Output directory for RGB image. Output directory will otherwise be the working directory",
		"Comma-separated list of RGB bands to use, e.g. -b 54,12,45",
		"Show this help",
		"Increase verbosity",
		"Do not generate new RGB image if an RGB image with the exact filename already exists",
		"Rescale intensities of image from mean - 2*stddev to mean + 2*stddev.",
		"Rescale intensitites of image to 0 to 1 using min and max intensities"
	};
	std::string usage_instructions = "\nUsage:\n" + std::string(argv[0]) + " [options] HYPERSPECTRAL_FILES\n";

	while (1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, short_options, long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'o': //output directory
				output_directory = std::string(optarg);
				break;
			case 'b': //comma-separated rgb bands
				rgb_bands_str = std::string(optarg);
				break;
			case 'h': //help
				getopt_long_show_help(usage_instructions.c_str(), long_options, short_options, option_descriptions);
				return 0;
				break;
			case OPT_IGNORE_EXISTING: //ignore existing RGB images
				ignore_existing_rgb_conversions = true;
				break;
			case 'v': //verbosity
				verbosity_level++;
				break;
			case OPT_RESCALE_WITH_STATISTICS:
				rescale_with_statistics = true;
				break;
			case OPT_RESCALE_TO_MAXMIN:
				rescale_to_maxmin = true;
				break;
		}
	}

	if (optind >= argc) {
		std::cerr << "Input filename missing, exiting." << std::endl;
		return -1;
	}

	//generate RGB conversions for each input hyperspectral file
	for (int i=optind; i < argc; i++) {
		std::string input_filename = std::string(argv[i]);
		if (verbosity_level > 0) {
			fprintf(stderr, "Processing %s...\n", input_filename.c_str());
		}

		//get basename of input filename
		char *input_filename_copy = strdup(input_filename.c_str());
		std::string output_filename(remove_file_ending(std::string(basename(input_filename_copy))) + ".png");
		free(input_filename_copy);

		//use specified output directory
		if (!output_directory.empty()) {
			output_filename = output_directory + "/" + output_filename;
		}

		//sanity check of whether output is the same as the input
		if (output_filename == input_filename) {
			fprintf(stderr, "Error in filename generation: Output filename the same as input filename.\n");
			continue;
		}

		//check whether we should ignore this file
		if (ignore_existing_rgb_conversions && (access(output_filename.c_str(), F_OK) != -1)) {
			if (verbosity_level > 0) {
				fprintf(stderr, "RGB conversion found, ignoring file.\n");
			}
			continue;
		}

		//read hyperspectral file
		struct hyspex_header header;
		hyperspectral_image_t hyperspectral_image;
		hyperspectral_err_t errcode = hyperspectral_read_header_and_image(input_filename.c_str(), &header, &hyperspectral_image);
		if (errcode != HYPERSPECTRAL_NO_ERR) {
			fprintf(stderr, "Error in reading %s: %s\n", input_filename.c_str(), hyperspectral_error_message(errcode));
			continue;
		}

		if (header.interleave != BIL_INTERLEAVE) {
			fprintf(stderr, "Unsupported interleave.\n");
			continue;
		}

		//get rgb bands from input
		std::vector<int> rgb_bands;
		if (!rgb_bands_str.empty()) {
			rgb_bands = split_string(rgb_bands_str, ",");
		} else {
			rgb_bands = header.default_bands;
		}

		std::vector<cv::Mat> channels = hyperspectral_to_rgb(&hyperspectral_image, rgb_bands, verbosity_level);
		if (channels.size() == 0) {
			continue;
		}

        double max_possible_value = powf(2, 16)-1;

		for (int i=0; i < channels.size(); i++) {
			//remove saturated pixels given that image is uint16_t (no effect on images already rescaled to [0,1]
			cv::Mat saturated_mask = channels[i] >= max_possible_value;
			channels[i].setTo(0, saturated_mask);

			//rescale image
			if (rescale_to_maxmin) {
				double min, max;
				cv::minMaxLoc(channels[i], &min, &max);
				channels[i] = (channels[i] - min)/(max - min);
			} else if (rescale_with_statistics) {
				cv::Scalar mean_out, stddev_out;
				cv::meanStdDev(channels[i], mean_out, stddev_out);

				float mean = mean_out[0];
				float stddev = stddev_out[0];

				float new_max = mean + 2*stddev;
				float new_min = mean - 2*stddev;

				channels[i].setTo(new_max, channels[i] > new_max);
				channels[i].setTo(new_min, channels[i] < new_min);
				channels[i] = (channels[i] - new_min)/(new_max - new_min);
			} else if (header.datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
                //no other rescaling suggested, raw data. Should rescale to 0, 1 from 0, 2^16 -1.
                channels[i] = channels[i]/max_possible_value;
            }

			channels[i] = channels[i]*255;
		}

		//reverse array since opencv expects BGR
		std::reverse(channels.begin(), channels.end());

		cv::Mat output_image;
		cv::merge(channels, output_image);
		cv::imwrite(output_filename, output_image);

		hyperspectral_free_data(&hyperspectral_image);
	}
}

#include <sys/mman.h>

std::vector<cv::Mat> hyperspectral_to_rgb(hyperspectral_image_t *image, std::vector<int> rgb_bands, int verbosity)
{
	struct hyspex_header header = image->header;
	int num_samples = header.samples;
	int num_lines = header.lines;
	int num_bands = header.bands;

	const int NUM_RGB_CHANNELS = 3;
	int num_channels = std::max<int>(rgb_bands.size(), NUM_RGB_CHANNELS);
	std::vector<cv::Mat> ret_channels;
	size_t element_size = hyperspectral_element_size(header); //size of each element in hyperspectral data array
	int cv_type;
	if (header.datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
		cv_type = CV_32FC1;
	} else if (header.datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
		cv_type = CV_16UC1;
	} else {
		fprintf(stderr, "Unsupported data type.\n");
		return std::vector<cv::Mat>();
	}

	char *channels = new char[image->header.samples*image->header.lines*num_channels*element_size]();
	char *image_data = hyperspectral_char(image);

	int pagesize = getpagesize();

	int k=0;
	std::vector<int> rgb_bands_sorted = rgb_bands;
	std::sort(rgb_bands_sorted.begin(), rgb_bands_sorted.end());

	//obtain data from file
	for (int i=0; i < header.lines; i++) {
		//prepare memory pattern for better file caching
		for (int j=0; j < num_channels; j++) {
			size_t address = i*header.samples*header.bands*element_size + rgb_bands_sorted[j]*header.samples*element_size;
			size_t length = element_size*header.samples;
			size_t aligned_address = (address/pagesize)*pagesize; //align address to pagesize
			length = length + (address - aligned_address);

			madvise(image_data + aligned_address, length, MADV_SEQUENTIAL | MADV_WILLNEED); //cache aggressively and sequentially
		}

		//read from file, in sorted band order
		for (int j=0; j < num_channels; j++) {
			if (verbosity > 0) print_progress_bar("Generating RGB image", k-1, k, (header.lines-1)*(num_channels-1));
			size_t dest_index = j*header.lines*header.samples + i*header.samples;
			size_t source_index = i*header.samples*header.bands + rgb_bands_sorted[j]*header.samples;
			memcpy(channels + dest_index*element_size, image_data + source_index*element_size, element_size*header.samples);
			k++;
		}
	}

	//convert to float type opencv arrays in correct RGB order
	for (int i=0; i < num_channels; i++) {
		int sorted_index = std::distance(rgb_bands_sorted.begin(), std::find(rgb_bands_sorted.begin(), rgb_bands_sorted.end(), rgb_bands[i]));

		cv::Mat channel_img;
		cv::Mat(num_lines, num_samples, cv_type, channels + sorted_index*header.samples*header.lines*element_size).convertTo(channel_img, CV_32FC1);
		ret_channels.push_back(channel_img);
	}
	return ret_channels;
}
