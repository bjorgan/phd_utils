#include <QApplication>
#include <string>
#include "getopt.h"
#include <sstream>
#include "imageViewer.h"
#include <hyperspectral/readimage.h>
#include <vector>
#include "string_helpers.h"
#include <iostream>
using namespace std;

void showHelp(){
	cerr << "Usage: hyview [OPTION]... [FILE]" << endl
		<< "Hyperspectral image viewer (BIL-interleaved ENVI image assumed)." << endl << endl
		<< "--help\t\t\t Show help" << endl << endl
		<< "Image subset arguments:" << endl
		<< "--startsample=START_PIXEL \t Start pixel (chosen pixel for showpixel)" << endl
		<< "--endsample=END_PIXEL \t End pixel" << endl
		<< "--startline=START_LINE \t Start line (chosen line for showpixel)" << endl
		<< "--endline=END_LINE \t End line" << endl
		<< "--startband=START_BAND \t Start band (chosen band for showpixel)" << endl
		<< "--endband=END_BAND \t End band" << endl
		<< "--minmax-rescaling \t Rescale to minmax instead of using mean/stdev" << endl;
}

void get_start_end(std::string input, int *ret_start, int *ret_end)
{
	std::vector<int> coordinates = split_string(input, ":");
	if (coordinates.size() != 2) {
		fprintf(stderr, "Insufficient number of arguments: %s\n", input.c_str());
		exit(1);
	}
	*ret_start = coordinates[0];
	*ret_end = coordinates[1];
}

enum opts{OPT_START_SAMPLE, OPT_END_SAMPLE, OPT_HELP, OPT_START_LINE, OPT_END_LINE, OPT_START_BAND, OPT_END_BAND, OPT_SAMPLES, OPT_LINES, OPT_BANDS, OPT_MINMAX_RESCALING};

int main(int argc, char *argv[]){
	//process options
	char shortopts[] = "";
	struct option longopts[] = {
		{"startsample",			required_argument,		0,	OPT_START_SAMPLE},
		{"endsample",			required_argument,		0,	OPT_END_SAMPLE},
		{"help",			no_argument,			0,	OPT_HELP},
		{"startline",			required_argument,		0,	OPT_START_LINE},
		{"endline",			required_argument,		0,	OPT_END_LINE},
		{"startband",			required_argument,		0,	OPT_START_BAND},
		{"endband",			required_argument,		0,	OPT_END_BAND},
		{"samples",			required_argument,		0,	OPT_SAMPLES},
		{"lines",			required_argument,		0,	OPT_LINES},
		{"bands",			required_argument,		0,	OPT_BANDS},
		{"minmax-rescaling",		no_argument,			0,	OPT_MINMAX_RESCALING},
		{0, 0, 0, 0}
	};

	int startline = 0;
	int endline = 0;
	int startsample = 0;
	int endsample = 0;
	int startband = 0;
	int endband = 0;

	int index;

	bool minmax_rescaling = false;

	while (true){
		int flag = getopt_long(argc, argv, shortopts, longopts, &index);
		switch (flag){
			case OPT_START_SAMPLE:
				startsample = strtod(optarg, NULL);
			break;

			case OPT_END_SAMPLE:
				endsample = strtod(optarg, NULL);
			break;

			case OPT_HELP: //help
				showHelp();
				exit(0);
			break;

			case OPT_START_LINE:
				startline = strtod(optarg, NULL);
			break;

			case OPT_END_LINE:
				endline = strtod(optarg, NULL);
			break;

			case OPT_START_BAND:
				startband = strtod(optarg, NULL);
			break;

			case OPT_END_BAND:
				endband = strtod(optarg, NULL);
			break;

			case OPT_SAMPLES:
				get_start_end(optarg, &startsample, &endsample);
			break;

			case OPT_LINES:
				get_start_end(optarg, &startline, &endline);
			break;

			case OPT_BANDS:
				get_start_end(optarg, &startband, &endband);
			break;

			case OPT_MINMAX_RESCALING:
				minmax_rescaling = true;
			break;
		}
		if (flag == -1){
			break;
		}
	}
	if (optind >= argc) {
		cerr << "Filename missing." << endl;
		exit(1);
	}
	std::string filename = std::string(argv[optind]);

	//read hyperspectral image header
	struct hyspex_header header;
	hyperspectral_read_header(filename.c_str(), &header);

	//configure image subsets
	if (!startline){
		startline = 0;
	}
	if (!endline){
		endline = header.lines;
	}
	if (!startsample){
		startsample = 0;
	}
	if (!endsample){
		endsample = header.samples;
	}
	if (!startband) {
		startband = 0;
	}
	if (!endband) {
		endband = header.bands;
	}

	int newLines = endline - startline;
	int newSamples = endsample - startsample;
	int newBands = endband - startband;

	struct image_subset subset = hyperspectral_generate_subset(header, startline, newLines, startsample, newSamples, startband, newBands);

	//read hyperspectral image
	float *data = hyperspectral_alloc_float(subset);
	hyperspectral_read_image(filename.c_str(), &header, subset, data);

	struct hyspex_header subset_header = hyperspectral_header_from_subset(header, subset);

	//start Qt app, display imageViewer widget
	QApplication app(argc, argv);
	ImageViewer viewer(subset, subset_header, data);

	if (minmax_rescaling) {
		viewer.minMaxRescaling();
	}

	viewer.show();

	return app.exec();
}
