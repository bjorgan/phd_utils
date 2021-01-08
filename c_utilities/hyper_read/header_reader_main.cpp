/**
 * CLI utility for printing hyperspectral image properties to stdout.
 **/

#include <string>
#include <iostream>
#include <getopt.h>
#include "getopt_long_helpers.h"
#include <hyperspectral/readimage.h>

int main(int argc, char *argv[])
{
	//getopt properties
	struct option long_options[] = {
		{"help",			no_argument,		0,	'h'},
		{"lines",			no_argument,		0,	'l'},
		{"samples",			no_argument,		0,	's'},
		{"bands",			no_argument,		0,	'b'},
		{"wavelengths",			no_argument,		0,	'w'},
		{0, 0, 0, 0}
	};
	char short_options[] = "hlsbw";
	const char *option_descriptions[] = {
		"Show help",
		NULL
	};
	std::string usage = "Usage: " + std::string(argv[0]) + " [OPTIONS] HYPERSPECTRAL_FILE\n";

	bool print_lines = false;
	bool print_samples = false;
	bool print_bands = false;
	bool print_wavelengths = false;

	int index;
	while (true){
		int flag = getopt_long(argc, argv, short_options, long_options, &index);
		switch (flag){
			case 'h': //help
				getopt_long_show_help(usage.c_str(), long_options, short_options, option_descriptions);
				exit(0);
			break;

			case 'l':
				print_lines = true;
			break;

			case 's':
				print_samples = true;
			break;

			case 'b':
				print_bands = true;
			break;

			case 'w':
				print_wavelengths = true;
			break;
		}
		if (flag == -1){
			break;
		}
	}

	if (optind >= argc) {
		fprintf(stderr, "Missing filename.\n");
		exit(1);
	}
    std::string filename = std::string(argv[optind]);

	struct hyspex_header header;
	hyperspectral_err_t ret = hyperspectral_read_header(filename.c_str(), &header);
	if (ret != HYPERSPECTRAL_NO_ERR){
		fprintf(stderr, "Error in reading image: %s, %s\n", filename.c_str(), hyperspectral_error_message(ret));
		exit(1);
	}

	if (print_lines) {
		std::cout << header.lines << " ";
	}

	if (print_samples) {
		std::cout << header.samples << " ";
	}

	if (print_bands) {
		std::cout << header.bands << " ";
	}

	if (print_wavelengths) {
		for (int i=0; i < header.bands; i++) {
			std::cout << header.wlens[i] << " ";
		}
	}
	std::cout << std::endl;
	return 0;
}

