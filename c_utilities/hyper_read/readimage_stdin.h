#ifndef READIMAGE_EXTENSIONS_H_DEFINED
#define READIMAGE_EXTENSIONS_H_DEFINED

#include <hyperspectral/readimage.h>

/**
 * Struct for controlling stdin as a hyperspectral stream.
 **/
typedef struct {
	///Whether stdin stream is open
	bool stream_closed;
	///Whether data should be saved to a larger array
	bool should_accumulate_data;
	///Hyperspectral header
	struct hyspex_header header;
	///Accumulated data
	char *accumulated_data;
	///Current line as last read from stdin
	char *curr_line_data;
	///Current line number in stdin (denotes next line to be read)
	int curr_line_in_stdin;
} hyperspectral_stdin_t;

/**
 * Options to stdin line reader.
 **/
enum hyperspectral_stdin_opt {
	STDIN_ACCUMULATE_DATA, //accumulate input data
	STDIN_NO_ACCUMULATION //don't accumulate input data
};

/**
 * Get line of data from stdin. Will spool until the specified line if current line
 * in stdin is less than the requested line.
 *
 * \param image stdin stream
 * \param line Line number. If corresponding to curr_line_in_stdin, this will read a new line of data from stdin. Otherwise, it will, if available, return the specified line from accumulated_data
 * \return NULL on failure, pointer to either curr_line_data or the appropriate line in accumulated_data on success
 **/
const char* hyperspectral_stdin_line_data(hyperspectral_stdin_t *image, int line);
hyperspectral_stdin_t *hyperspectral_stdin_create(struct hyspex_header header, int hyperspectral_stdin_options = STDIN_NO_ACCUMULATION);

/**
 * General hyperspectal reader, either from file or stdin.
 **/
typedef struct {
	///Hyperspectral header, read from somewhere else
	struct hyspex_header header;
	///Whether we should read from stdin
	bool read_from_stdin;
	///stdin reader
	hyperspectral_stdin_t *stdin_stream;
	///File reader
	hyperspectral_image_t *file_stream;
} hyperspectral_reader_t;

/**
 * Create general hyperspectral reader.
 *
 * \param header Hyperspectral header
 * \param filename Filename. If empty, it is assumed that we will read from stdin
 * \param stdin_options Input options to stdin reading
 **/
hyperspectral_reader_t *hyperspectral_create_reader(struct hyspex_header header, std::string filename, int stdin_options = STDIN_NO_ACCUMULATION);

/**
 * Get hyperspectral line of data from hyperspectral reader. Can jump back and forth if we are reading
 * from file or reading from stdin and accumulating data, will otherwise have to do it incrementally, line by line, increasing the line number on every read
 *
 * \param reader Hyperspectral reader
 * \param line Line in image
 * \return Pointer to line on success, NULL on failure
 **/
const char *hyperspectral_line_data(hyperspectral_reader_t *reader, int line);

/**
 * Convert hyperspectral raw array to float array.
 *
 * \param header Hyperspectral header
 * \param line Input raw line
 * \param output Output data
 **/
void hyperspectral_line_to_float(struct hyspex_header header, const char *line, float *output);

/**
 * Read line from hyperspectral reader and copy and convert data to float array.
 *
 * \param reader Hyperspectral reader
 * \param line Line number
 * \param output Output data array
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_line_data_float(hyperspectral_reader_t *reader, int line, float *output);

#endif
