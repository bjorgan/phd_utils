#include "readimage.h"
#include <boost/regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <numeric>
#include <sys/mman.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libexplain/mmap.h>

///Max size for c strings
const int MAX_NUM_CHARS = 512;
///Maximum assumed number of characters in header file
const int MAX_HEADER_FILE_SIZE = 4000;

/**
 * For defining whether to look for a scalar or vector type property value in the header.
 **/
enum header_property_type {
	///Scalar value
	SCALAR_VALUE,
	///Vector value
	ARRAY_VALUE
};

/**
 * Look for property in header text and return property value as string, given format `property = value` (with arbitrary whitespaces).
 *
 * \param hdrText Full header text
 * \param property Property to look for
 * \param header_property_type Property type, scalar or vector
 * \return Value associated with property
 **/
std::string header_extract_property(std::string hdrText, std::string property, enum header_property_type header_property_type = SCALAR_VALUE);

/**
 * Look for array property in header text and return array values as a single string.
 **/
std::string header_extract_array_property(std::string hdrText, std::string property);

/**
 * Split wavelength string (`{wlen1, wlen2, ...}` (or whitespace-separated)) into wavelength list.
 *
 * \param bands Number of wavelengths
 * \param wavelengthStr Wavelenght string
 * \return List of wavelengths
 **/
std::vector<float> split_wavelengths(int bands, std::string wavelengthStr);

/**
 * Remove file ending from filename.
 *
 * \param filename Filename
 * \return Filename with removed file ending (assumed to be last .)
 **/
std::string remove_file_ending(std::string filename);

/**
 * Read full hyperspectral file into a float array using conventional fread.
 *
 * \param filename Filename
 * \param header Header
 * \param data Output data
 **/
hyperspectral_err_t hyperspectral_read_image_sequential(const char *filename, struct hyspex_header *header, float *data);

/**
 * Read hyperspectral file into a float array using mmap and subsetting.
 *
 * \param filename Filename
 * \param header Header
 * \param subset Image subset
 * \param data Output data
 **/
hyperspectral_err_t hyperspectral_read_image_mmap(const char *filename, struct hyspex_header *header, struct image_subset subset, float *data);

size_t hyperspectral_element_size(struct hyspex_header header)
{
	int datatype = header.datatype;
	size_t element_bytes = 0;
	if (datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
		element_bytes = sizeof(float);
	} else if (datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
		element_bytes = sizeof(uint16_t);
	}
	return element_bytes;
}

#include <algorithm>

hyperspectral_err_t hyperspectral_header_from_string(std::string header_text, struct hyspex_header *header)
{
	//lowercasify full header text
	std::transform(header_text.begin(), header_text.end(), header_text.begin(), ::tolower);

	//reset arrays
	header->wlens.clear();
	header->default_bands.clear();

	//extract properties from header file text
	std::string samples = header_extract_property(header_text, "samples");
	std::string bands = header_extract_property(header_text, "bands");
	std::string lines = header_extract_property(header_text, "lines");
	std::string wavelengths = header_extract_property(header_text, "wavelength", ARRAY_VALUE);
	std::string hdrOffset = header_extract_property(header_text, "header offset");
	std::string interleave = header_extract_property(header_text, "interleave");
	std::string datatype = header_extract_property(header_text, "data type");
	std::string default_bands_str = header_extract_property(header_text, "default bands", ARRAY_VALUE);

	if ((samples.empty()) || (bands.empty()) || (lines.empty()) || (hdrOffset.empty()) || (interleave.empty()) || (datatype.empty())) {
		return HYPERSPECTRAL_HDR_PROPERTY_NOT_FOUND;
	}

	if (default_bands_str.empty()) {
		default_bands_str = "{0,0,0}";
	}

	//convert strings to values
	header->bands = strtod(bands.c_str(), NULL);
	header->lines = strtod(lines.c_str(), NULL);
	header->samples = strtod(samples.c_str(), NULL);
	header->offset = strtod(hdrOffset.c_str(), NULL);
	header->wlens = split_wavelengths(header->bands, wavelengths);
	header->datatype = strtod(datatype.c_str(), NULL);

	if (hyperspectral_element_size(*header) == 0) {
		return HYPERSPECTRAL_DATATYPE_UNSUPPORTED;
	}

	std::vector<float> default_bands = split_wavelengths(3, default_bands_str);
	for (int i=0; i < default_bands.size(); i++) {
		header->default_bands.push_back(default_bands[i]);
	}

	if (interleave == "bil") {
		header->interleave = BIL_INTERLEAVE;
	} else {
		return HYPERSPECTRAL_INTERLEAVE_UNSUPPORTED;
	}

	return HYPERSPECTRAL_NO_ERR;
}

hyperspectral_err_t hyperspectral_read_header_from_binary_stream(FILE *stream, struct hyspex_header *header)
{
	size_t num_bytes_in_header;
	int read_bytes = fread(&num_bytes_in_header, sizeof(size_t), 1, stream);
	if (read_bytes <= 0) {
		return HYPERSPECTRAL_FILE_READING_ERROR;
	}

	char *header_text;
	try {
		header_text = new char[num_bytes_in_header+1]();
	} catch (const std::exception &ex) {
		return HYPERSPECTRAL_FILE_READING_ERROR;
	}
	read_bytes = fread(header_text, sizeof(char), num_bytes_in_header, stream);
	if (read_bytes <= 0) {
		return HYPERSPECTRAL_FILE_READING_ERROR;
	}

	header_text[num_bytes_in_header] = '\0';

	std::string header_text_str = std::string(header_text);
	delete [] header_text;

	return hyperspectral_header_from_string(header_text_str, header);
}

hyperspectral_err_t hyperspectral_read_header_from_stdin(struct hyspex_header *header)
{
	hyperspectral_err_t errcode = hyperspectral_read_header_from_binary_stream(stdin, header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		return HYPERSPECTRAL_STDIN_FAILED;
	}
}

#include <sstream>

hyperspectral_err_t hyperspectral_read_header(const char* filename, struct hyspex_header *header)
{
	if (strlen(filename) == 0) {
		return hyperspectral_read_header_from_stdin(header);
	}

	//find base filename
	std::string baseName = remove_file_ending(std::string(filename));

	//open and read header file
	//append .hdr to filename
	std::string hdrName = baseName + ".hdr";

	FILE *fp = fopen(hdrName.c_str(), "rt");
	if (fp == NULL) {
		return HYPERSPECTRAL_HEADER_FILE_NOT_FOUND;
	}

	//read header text from file
	std::ifstream in_file;
	in_file.open(hdrName);

	std::stringstream str_stream;
	str_stream << in_file.rdbuf();

	std::string hdrText = str_stream.str();

	return hyperspectral_header_from_string(hdrText, header);
}

hyperspectral_err_t hyperspectral_read_image_mmap(const char *filename, struct hyspex_header *header, struct image_subset subset, float *data)
{
	hyperspectral_image_t hyperspectral_map;
	hyperspectral_err_t err = hyperspectral_read_image(filename, *header, &hyperspectral_map);
	if (err != HYPERSPECTRAL_NO_ERR) {
		return err;
	}

	//skip header and lines we do not want
	size_t elementBytes = hyperspectral_element_size(*header);
	size_t skipBytes = subset.start_line*header->samples*header->bands*elementBytes + header->offset;
	int num_samples = 0, num_lines = 0, num_bands = 0;
	hyperspectral_get_size(subset, &num_lines, &num_bands, &num_samples);

	//read in line by line from mmap
	for (int i=0; i < num_lines; i++) {
		char *line = (char*)hyperspectral_map.mmapped_file + skipBytes + i*elementBytes*header->samples*header->bands;
		for (int k=0; k < num_bands; k++) {
			for (int j=subset.start_sample; j < subset.end_sample; j++) {
				float val;
				int position = (subset.start_band + k)*header->samples + j;
				if (header->datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
					val = ((float*)line)[position];
				} else if (header->datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
					val = (((uint16_t*)line)[position])*1.0f;
				}

				data[i*num_bands*num_samples + k*num_samples + j-subset.start_sample] = val;
			}
		}
	}


	hyperspectral_free_data(&hyperspectral_map);
	return HYPERSPECTRAL_NO_ERR;
}

hyperspectral_err_t hyperspectral_read_image_sequential(const char *filename, struct hyspex_header *header, float *data)
{
	FILE *fp = fopen(filename, "rb");
	if (fp == NULL){
		return HYPERSPECTRAL_IMAGE_FILE_NOT_FOUND;
        }
	fseek(fp, header->offset, SEEK_SET);

        //read in line by line directly using fread
	size_t elementBytes = hyperspectral_element_size(*header);
	int num_pixels_in_line = header->bands*header->samples;
	char *in_line = new char[elementBytes*num_pixels_in_line]();
	for (int i=0; i < header->lines; i++) {
		int size_read = fread(in_line, elementBytes, num_pixels_in_line, fp);
		if (size_read == 0) {
			delete [] in_line;
			fclose(fp);
			return HYPERSPECTRAL_FILE_READING_ERROR;
		}

		float *out_line = data + i*num_pixels_in_line;
		if (header->datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
			memcpy(out_line, in_line, num_pixels_in_line*elementBytes);;
		} else {
			for (int j=0; j < num_pixels_in_line; j++) {
				float val = 0.0f;
				if (header->datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
					val = (((uint16_t*)in_line)[j])*1.0f;
				}
				out_line[j] = val;
			}
		}
	}
	delete [] in_line;
	fclose(fp);
}

hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header *header, struct image_subset subset, float *data)
{
	return hyperspectral_read_image_mmap(filename, header, subset, data);
}

hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header *header, float *data)
{
	return hyperspectral_read_image_sequential(filename, header, data);
}

float *hyperspectral_float_line(hyperspectral_image_t *data, int line)
{
	return (float*)hyperspectral_char_line(data, line);
}

float *hyperspectral_float(hyperspectral_image_t *data)
{
	return (float*)hyperspectral_char(data);
}

uint16_t *hyperspectral_uint16_line(hyperspectral_image_t *data, int line)
{
	return (uint16_t*)hyperspectral_char_line(data, line);
}

uint16_t *hyperspectral_uint16(hyperspectral_image_t *data)
{
	return (uint16_t*)hyperspectral_char(data);
}

char *hyperspectral_char(hyperspectral_image_t *data)
{
	return (char*)data->mmapped_file + data->header.offset;
}

char *hyperspectral_char_line(hyperspectral_image_t *data, int line)
{
	return (char*)data->mmapped_file + data->header.offset + line*data->header.samples*data->header.bands*hyperspectral_element_size(data->header);
}

void hyperspectral_free_data(hyperspectral_image_t *data)
{
	size_t elementBytes = hyperspectral_element_size(data->header);
	size_t num_bytes = elementBytes*data->header.lines*data->header.bands*data->header.samples;
	munmap(data->mmapped_file, data->mapped_size);
	close(data->file_descriptor);
}

hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header header, hyperspectral_image_t *data)
{
	data->header = header;
	size_t elementBytes = hyperspectral_element_size(header);

	data->file_descriptor = open(filename, O_RDONLY);
	if (data->file_descriptor < 0) {
		return HYPERSPECTRAL_IMAGE_FILE_NOT_FOUND;
	}

	//get size of file, compare against expected size
	size_t num_bytes = elementBytes*header.lines*header.bands*header.samples;

	struct stat stbuf;
	if ((fstat(data->file_descriptor, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode))) {
		return HYPERSPECTRAL_MMAP_FAILED;
	}
	off_t file_size = stbuf.st_size;
	if (file_size < num_bytes + header.offset) {
		fprintf(stderr, "Physical file size less than expected file size.\n");
		return HYPERSPECTRAL_MMAP_FAILED;
	}

	//get mmapping of file
	data->mapped_size = num_bytes + header.offset;
	data->mmapped_file = mmap(NULL, data->mapped_size, PROT_READ, MAP_PRIVATE, data->file_descriptor, 0);
	if (data->mmapped_file == MAP_FAILED) {
		fprintf(stderr, "mmap failed with error message %s\n", explain_mmap(NULL, data->mapped_size, PROT_READ, MAP_PRIVATE, data->file_descriptor, 0));
		return HYPERSPECTRAL_MMAP_FAILED;
	}

	return HYPERSPECTRAL_NO_ERR;
}

hyperspectral_err_t hyperspectral_read_header_and_image(const char *filename, struct hyspex_header *header, hyperspectral_image_t *data)
{
	hyperspectral_err_t errcode = hyperspectral_read_header(filename, header);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		return errcode;
	}

	errcode = hyperspectral_read_image(filename, *header, data);
	if (errcode != HYPERSPECTRAL_NO_ERR) {
		return errcode;
	}

	return HYPERSPECTRAL_NO_ERR;
}


std::string getMatch(std::string input_string, regmatch_t *matchArray, int matchNum)
{
	int start = matchArray[matchNum].rm_so;
	int end = matchArray[matchNum].rm_eo;
	return input_string.substr(start, end - start);
}

std::string header_extract_property(std::string hdrText, std::string property)
{
	return header_extract_property(hdrText, property, SCALAR_VALUE);
}

std::string header_extract_array_property(std::string hdrText, std::string property)
{
	return header_extract_property(hdrText, property, ARRAY_VALUE);
}

std::string header_extract_property(std::string hdrText, std::string property, enum header_property_type header_property_type)
{
	//find the input property using regex
	regex_t propertyMatch;
	int numMatch = 2;
	regmatch_t *matchArray = (regmatch_t*)malloc(sizeof(regmatch_t)*numMatch);

	char regexExpr[MAX_NUM_CHARS] = "";
	strcat(regexExpr, property.c_str());

	if (header_property_type == SCALAR_VALUE) {
		//match full value
		strcat(regexExpr, "\\s*=\\s*([0-9|,| |.|a-z]+)"); //property followed by = and a set of number, commas or spaces
	} else {
		//match only on start of array due to multiline difficulty
		strcat(regexExpr, "\\s*=\\s*{"); //property followed by = {and a set of number, commas or spaces, newlines}
	}
	int retcode = regcomp(&propertyMatch, regexExpr, REG_EXTENDED | REG_NEWLINE | REG_PERL);
	int match = regexec(&propertyMatch, hdrText.c_str(), numMatch, matchArray, 0);
	if (match != 0) {
		regfree(&propertyMatch);
		free(matchArray);
		return std::string();
	}

	std::string retVal;
	if (header_property_type == SCALAR_VALUE) {
		//get full scalar value from regex match
		retVal = getMatch(hdrText, matchArray, 1);
	} else {
		//get start posision of array, and extract until first }
		int array_start = matchArray[0].rm_eo;
		int array_end = 0;
		std::string substring = hdrText.substr(array_start, hdrText.length() - array_start);
		for (int i=0; i < substring.length(); i++) {
			if (substring[i] == '}') {
				array_end = i;
				break;
			}
		}

		retVal = substring.substr(0, array_end);
	}

	//cleanup
	regfree(&propertyMatch);
	free(matchArray);

	return retVal;
}

std::string remove_file_ending(std::string filename)
{
	regex_t filenameMatch;
	int numMatch = 2;
	regmatch_t *matchArray = (regmatch_t*)malloc(sizeof(regmatch_t)*numMatch);
	int retcode = regcomp(&filenameMatch, "(.*)[.].*$", REG_EXTENDED);
	int match = regexec(&filenameMatch, filename.c_str(), numMatch, matchArray, 0);

	std::string baseName = filename;
	if (match == 0) {
		baseName = getMatch(filename, matchArray, 1);
	}

	regfree(&filenameMatch);
	free(matchArray);
	return baseName;
}

std::vector<float> split_wavelengths(int bands, std::string wavelengthStr)
{
	std::vector<float> retWlens;
	bool useStandardValues = false;
	if (!wavelengthStr.empty()) {
		//prepare regex
		regex_t numberMatch;
		int numMatch = 2;
		regmatch_t *matchArray = (regmatch_t*)malloc(sizeof(regmatch_t)*numMatch);
		char regexExpr[MAX_NUM_CHARS] = "([0-9|.]+)[,| |}]*";
		int retcode = regcomp(&numberMatch, regexExpr, REG_EXTENDED);

		//find start of number sequence
		int currStart = wavelengthStr.find_first_of('{') + 1;
		std::string substr = wavelengthStr.substr(currStart);

		//go through all bands
		for (int i=0; i < bands; i++) {
			substr = wavelengthStr.substr(currStart);

			//extract wavelength
			if (regexec(&numberMatch, substr.c_str(), numMatch, matchArray, 0)) {
				useStandardValues = true;
				break;
			}

			std::string match = getMatch(substr.c_str(), matchArray, 1);
			retWlens.push_back(strtod(match.c_str(), NULL));

			//move to next
			currStart = currStart + matchArray[0].rm_eo;
		}
		regfree(&numberMatch);
		free(matchArray);
	} else {
		useStandardValues = true;
	}

	if (useStandardValues) {
		retWlens.clear();
		for (int i=0; i < bands; i++) {
			retWlens.push_back(i);
		}
	}
	return retWlens;
}

const char *hyperspectral_error_message(hyperspectral_err_t error_code)
{
	switch (error_code) {
		case HYPERSPECTRAL_NO_ERR:
			return "No error.";
		case HYPERSPECTRAL_HEADER_FILE_NOT_FOUND:
			return "Header file not found.";
		case HYPERSPECTRAL_IMAGE_FILE_NOT_FOUND:
			return "Image file not found.";
		case HYPERSPECTRAL_HDR_PROPERTY_NOT_FOUND:
			return "Could not find a required header file property.";
		case HYPERSPECTRAL_INTERLEAVE_UNSUPPORTED:
			return "Interleave not supported.";
		case HYPERSPECTRAL_DATATYPE_UNSUPPORTED:
			return "Datatype not supported.";
		case HYPERSPECTRAL_FILE_READING_ERROR:
			return "General hyperspectral file reading error.";
		case HYPERSPECTRAL_MMAP_FAILED:
			return "mmap() failed to map hyperspectral image to array.";
		case HYPERSPECTRAL_STDIN_FAILED:
			return "Failed to read hyperspectral information from stdin.";
	}
}


#include <fstream>
#include <iostream>
#include <math.h>
#include <cstring>
#include <sstream>
using namespace std;

struct image_subset hyperspectral_generate_subset(struct hyspex_header header, int start_line, int num_lines, int start_sample, int num_samples, int start_band, int num_bands)
{
	struct image_subset subset;
	subset.start_line = start_line;
	subset.end_line = start_line + num_lines;
	if (subset.end_line > header.lines) subset.end_line = header.lines;
	if (subset.start_line < 0) subset.start_line = 0;
	if (subset.start_line > subset.end_line) subset.start_line = subset.end_line;

	subset.start_sample = start_sample;
	subset.end_sample = start_sample + num_samples;
	if (subset.end_sample > header.samples) subset.end_sample = header.samples;
	if (subset.start_sample < 0) subset.start_sample = 0;
	if (subset.start_sample > subset.end_sample) subset.start_sample = subset.end_sample;

	subset.start_band = start_band;
	subset.end_band = start_band + num_bands;
	if (subset.end_band > header.bands) subset.end_band = header.bands;
	if (subset.start_band < 0) subset.start_band = 0;
	if (subset.start_band > subset.end_band) subset.start_band = subset.end_band;
	return subset;
}

float *hyperspectral_alloc_float(struct hyspex_header header)
{
	return new float[header.lines*header.bands*header.samples]();
}

void hyperspectral_get_size(struct image_subset image_subset, int *num_lines, int *num_bands, int *num_samples)
{
	*num_bands = image_subset.end_band - image_subset.start_band;
	*num_lines = image_subset.end_line - image_subset.start_line;
	*num_samples = image_subset.end_sample - image_subset.start_sample;
}

float *hyperspectral_alloc_float(struct image_subset image_subset)
{
	int num_samples = 0, num_lines = 0, num_bands = 0;
	hyperspectral_get_size(image_subset, &num_lines, &num_bands, &num_samples);
	return new float[num_lines*num_samples*num_bands]();
}

struct image_subset hyperspectral_generate_subset(hyspex_header header)
{
	struct image_subset subset = {0};
	subset.end_line = header.lines;
	subset.end_band = header.bands;
	subset.end_sample = header.samples;
	return subset;
}

const char* header_interleave_string(struct hyspex_header header)
{
	switch (header.interleave) {
		case BIL_INTERLEAVE:
			return "bil";
		case BIP_INTERLEAVE:
			return "bip";
		case BSQ_INTERLEAVE:
			return "bsq";
	}
}

std::string hyperspectral_header_to_string(struct hyspex_header header)
{
	std::ostringstream header_output;
	header_output << "ENVI" << endl;
	header_output << "samples = " << header.samples << endl;
	header_output << "lines = " << header.lines << endl;
	header_output << "bands = " << header.bands << endl;
	header_output << "header offset = " << header.offset << endl;
	header_output << "file type = ENVI Standard" << endl;
	header_output << "data type = " << header.datatype << endl;
	header_output << "interleave = " << header_interleave_string(header) << endl;
	header_output << "default bands = {";
	for (int i=0; i < header.default_bands.size(); i++) {
		header_output << header.default_bands[i];
		if (i < header.default_bands.size()-1) {
			header_output << ",";
		}
	}
	header_output << "}" << endl;
	header_output << "byte order = 0" << endl;
	header_output << "wavelength = {";
	for (int i=0; i < header.wlens.size(); i++) {
		header_output << header.wlens[i] << " ";
	}
	header_output << "}" << endl;
	return header_output.str();
}

struct hyspex_header hyperspectral_header_from_sizes(int bands, int samples, int lines, std::vector<float> wlens)
{
	struct hyspex_header header;
	header.samples = samples;
	header.bands = bands;
	header.lines = lines;
	header.wlens = wlens;
	header.default_bands.push_back(55);
	header.default_bands.push_back(41);
	header.default_bands.push_back(12);
	header.interleave = BIL_INTERLEAVE;
	header.datatype = HYPERSPECTRAL_DATATYPE_FLOAT;
	header.offset = 0;
	return header;
}

std::string hyperspectral_header_to_string(int bands, int samples, int lines, std::vector<float> wlens)
{
	return hyperspectral_header_to_string(hyperspectral_header_from_sizes(bands, samples, lines, wlens));
}

/**
 * Used for specifying whether the header should be written at the start of the file stream or not.
 **/
enum hyperspectral_file_option {WRITE_HEADER_TO_IMAGE_STREAM, DO_NOT_WRITE_HEADER_TO_IMAGE_STREAM};

/**
 * Generate hyperspectral_file_t from a FILE* stream (already opened). Used only internally, don't expose
 * to the public API. Used for having a common function for stdout and actual file streams, and
 * for proper testing in the testing framework.
 *
 * \param stream File stream
 * \param header Hyperspectral header
 * \param option Whether header struct should be written to the file stream. NB! FIXME! The final file header, if
 * written to file, will not contain an offset, so if this is an actual image, the image will be broken.
 * \return hyperspectral_file_t instance
 **/
hyperspectral_file_t hyperspectral_file_from_stream(FILE *stream, struct hyspex_header header, hyperspectral_file_option option = DO_NOT_WRITE_HEADER_TO_IMAGE_STREAM)
{
	hyperspectral_file_t ret_file;
	ret_file.filename = std::string();
	ret_file.header = header;
	ret_file.written_bytes = 0;
	ret_file.image_file = stream;

	if (option == WRITE_HEADER_TO_IMAGE_STREAM) {
		std::string header_string = hyperspectral_header_to_string(header);
		size_t num_chars_in_header = strlen(header_string.c_str());
		hyperspectral_write_to_file(&ret_file, sizeof(size_t), (char*)&num_chars_in_header);
		hyperspectral_write_to_file(&ret_file, num_chars_in_header, header_string.c_str());
	}

	return ret_file;
}

/**
 * Generate hyperspectral_file_t from a FILE* stream (already opened). Used only internally, don't expose
 * to the public API. Used for having a common function for stdout and actual file streams, and
 * for proper testing in the testing framework.
 *
 * \param stream File stream
 * \param lines Number of lines (used in header when using option WRITE_HEADER_TO_IMAGE_STREAM)
 * \param bands Number of bands
 * \param samples Number of bands
 * \param wlens Wavelengths
 * \param option Whether header struct should be written to the file stream. NB! FIXME! The final file header, if
 * written to file, will not contain an offset, so if this is an actual image, the image will be broken.
 * \return hyperspectral_file_t instance
 **/
hyperspectral_file_t hyperspectral_file_from_stream(FILE *stream, int lines, int bands, int samples, std::vector<float> wlens, hyperspectral_file_option option = DO_NOT_WRITE_HEADER_TO_IMAGE_STREAM)
{
	struct hyspex_header header = hyperspectral_header_from_sizes(bands, samples, lines, wlens);
	return hyperspectral_file_from_stream(stream, header, option);
}

const int NUM_LINES_UNKNOWN = -1;

hyperspectral_file_t hyperspectral_open_write_file(std::string filename, struct hyspex_header header)
{
	std::string out_filename = std::string(filename) + IMAGE_EXTENSION;
	FILE *stream = fopen(out_filename.c_str(), "wb");

	hyperspectral_file_t ret_file = hyperspectral_file_from_stream(stream, header);
	ret_file.filename = filename;
	return ret_file;
}

hyperspectral_file_t hyperspectral_open_write_file(std::string filename, int bands, int samples, std::vector<float> wlens)
{
	struct hyspex_header header = hyperspectral_header_from_sizes(bands, samples, NUM_LINES_UNKNOWN, wlens);
	return hyperspectral_open_write_file(filename, header);
}

hyperspectral_file_t hyperspectral_open_stdout_file(struct hyspex_header header)
{
	return hyperspectral_file_from_stream(stdout, header, WRITE_HEADER_TO_IMAGE_STREAM);
}

hyperspectral_file_t hyperspectral_open_stdout_file(int lines, int bands, int samples, std::vector<float> wlens)
{
	struct hyspex_header header = hyperspectral_header_from_sizes(bands, samples, lines, wlens);
	return hyperspectral_open_stdout_file(header);
}


size_t hyperspectral_write_to_file(hyperspectral_file_t *file, size_t num_bytes, const char *data)
{
	size_t num_written_bytes = fwrite(data, sizeof(char), num_bytes, file->image_file);
	file->written_bytes += num_bytes;
	return num_written_bytes;
}

void hyperspectral_free(hyperspectral_file_t *file)
{
	if (!file->filename.empty()) {
		fclose(file->image_file);
	}
}

void hyperspectral_close_write_file(hyperspectral_file_t *file)
{
	hyperspectral_free(file);

	if (!file->filename.empty()) {
		file->header.lines = file->written_bytes/(file->header.samples*file->header.bands*hyperspectral_element_size(file->header));
		hyperspectral_write_header(file->filename.c_str(), file->header);
	}
}

void hyperspectral_write_header(const char *filename, int bands, int samples, int lines, std::vector<float> wlens)
{
	std::string out_filename = std::string(filename) + HEADER_EXTENSION;

	std::ofstream header_file(out_filename.c_str());
	header_file << hyperspectral_header_to_string(bands, samples, lines, wlens);
	header_file.close();
}

void hyperspectral_write_image(const char *filename, int bands, int samples, int lines, float *data)
{
	hyperspectral_file_t file = hyperspectral_open_write_file(filename, bands, samples, std::vector<float>());
	hyperspectral_write_to_file(&file, lines*bands*samples*sizeof(float), (char*)data);
	hyperspectral_free(&file);
}

hyperspectral_file_t hyperspectral_open(std::string filename, int lines, int bands, int samples, std::vector<float> wlens)
{
	if (filename.empty()) {
		return hyperspectral_open_stdout_file(lines, bands, samples, wlens);
	} else {
		return hyperspectral_open_write_file(filename, bands, samples, wlens);
	}
}

hyperspectral_file_t hyperspectral_open(std::string filename, struct hyspex_header header)
{
	if (filename.empty()) {
		return hyperspectral_open_stdout_file(header);
	} else {
		return hyperspectral_open_write_file(filename, header);
	}
}

void hyperspectral_write_image(std::string filename, struct hyspex_header header, char *data)
{
	hyperspectral_file_t file = hyperspectral_open_write_file(filename, header.bands, header.samples, header.wlens);
	hyperspectral_write_to_file(&file, hyperspectral_element_size(header)*header.samples*header.lines*header.bands, data);
	hyperspectral_close_write_file(&file);
}

void hyperspectral_write_header(std::string filename, struct hyspex_header header)
{
	std::string out_filename = filename + HEADER_EXTENSION;
	std::ofstream header_file(out_filename.c_str());
	header_file << hyperspectral_header_to_string(header);
	header_file.close();
}

struct hyspex_header hyperspectral_header_from_subset(struct hyspex_header header, struct image_subset subset)
{
	struct hyspex_header ret_header = header;
	int samples, lines, bands;
	hyperspectral_get_size(subset, &lines, &bands, &samples);
	ret_header.samples = samples;
	ret_header.lines = lines;
	ret_header.bands = bands;
	vector<float>::const_iterator first = header.wlens.begin() + subset.start_band;
	vector<float>::const_iterator last = header.wlens.begin() + subset.end_band;
	std::vector<float> reduced_wlens(first, last);
	ret_header.wlens = reduced_wlens;
	return ret_header;
}
