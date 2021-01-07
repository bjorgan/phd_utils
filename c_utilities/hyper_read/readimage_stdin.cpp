#include "readimage_stdin.h"
#include <cstring>

const char* hyperspectral_stdin_line_data(hyperspectral_stdin_t *image, int line)
{
	size_t line_size = image->header.bands*image->header.samples;
	size_t line_bytes = line_size*hyperspectral_element_size(image->header);

	if ((line == image->curr_line_in_stdin) && !image->stream_closed) {
		//read from stdin
		int num_elements = fread(image->curr_line_data, hyperspectral_element_size(image->header), line_size, stdin);
		if (num_elements < line_size) {
			fprintf(stderr, "Error in stdin reading.\n");
			exit(1);
		}

		if (image->should_accumulate_data) {
			memcpy(image->accumulated_data + image->curr_line_in_stdin*line_bytes, image->curr_line_data, line_bytes);
		}

		image->curr_line_in_stdin++;
		if (image->curr_line_in_stdin >= image->header.lines) {
			image->stream_closed = true;
		}
		return image->curr_line_data;
	} else if ((line < image->curr_line_in_stdin) && image->should_accumulate_data) {
		//look up in existing data (if any)
		return image->accumulated_data + line*line_bytes;
	} else if ((line > image->curr_line_in_stdin) && !image->stream_closed) {
		//seek until desired line
		for (int i=image->curr_line_in_stdin; i < line; i++) {
			if (hyperspectral_stdin_line_data(image, i) == NULL) {
				return NULL;
			}
		}
		return hyperspectral_stdin_line_data(image, line);
	} else {
		return NULL;
	}
}

hyperspectral_stdin_t *hyperspectral_stdin_create(struct hyspex_header header, int hyperspectral_stdin_options)
{
	hyperspectral_stdin_t *image = new hyperspectral_stdin_t;
	image->stream_closed = false;
	image->should_accumulate_data = (hyperspectral_stdin_options == STDIN_ACCUMULATE_DATA);
	image->header = header;
	image->curr_line_data = new char[hyperspectral_element_size(header)*header.samples*header.bands]();
	if (image->should_accumulate_data) {
		image->accumulated_data = new char[hyperspectral_element_size(header)*header.lines*header.samples*header.bands]();
	} else {
		image->accumulated_data = NULL;
	}
	image->curr_line_in_stdin = 0;
	return image;
}

hyperspectral_reader_t *hyperspectral_create_reader(struct hyspex_header header, std::string filename, int stdin_options)
{
	hyperspectral_reader_t *reader = new hyperspectral_reader_t;
	reader->header = header;
	if (filename.empty()) {
		reader->read_from_stdin = true;
		reader->stdin_stream = hyperspectral_stdin_create(header, stdin_options);
	} else {
		reader->read_from_stdin = false;
		reader->file_stream = new hyperspectral_image_t;
		hyperspectral_err_t errcode = hyperspectral_read_image(filename.c_str(), header, reader->file_stream);
		if (errcode != HYPERSPECTRAL_NO_ERR) {
			fprintf(stderr, "Error opening %s: %s\n", filename.c_str(), hyperspectral_error_message(errcode));
		}
	}
	return reader;
}

const char *hyperspectral_line_data(hyperspectral_reader_t *reader, int line)
{
	const char *data = NULL;
	if (reader->read_from_stdin) {
		data = hyperspectral_stdin_line_data(reader->stdin_stream, line);
	} else {
		data = hyperspectral_char_line(reader->file_stream, line);
	}
}

void hyperspectral_line_to_float(struct hyspex_header header, const char *line, float *output)
{
	if (header.datatype == HYPERSPECTRAL_DATATYPE_FLOAT) {
		memcpy(output, line, sizeof(float)*header.samples*header.bands);
	} else if (header.datatype == HYPERSPECTRAL_DATATYPE_UINT16_T) {
		for (int i=0; i < header.samples*header.bands; i++) {
			output[i] = ((uint16_t*)line)[i];
		}
	}
}

hyperspectral_err_t hyperspectral_line_data_float(hyperspectral_reader_t *reader, int line, float *output)
{
	const char *line_data = hyperspectral_line_data(reader, line);
	if (line_data == NULL) {
		return HYPERSPECTRAL_FILE_READING_ERROR;
	}
	hyperspectral_line_to_float(reader->header, line_data, output);
	return HYPERSPECTRAL_NO_ERR;
}
