#ifndef READIMAGE_H_DEFINED
#define READIMAGE_H_DEFINED
#include <vector>
#include <string>
#include <stdint.h>

const std::string IMAGE_EXTENSION = ".img";
const std::string HEADER_EXTENSION = ".hdr";

/**
 * Errors.
 **/
enum hyperspectral_err_t {
	///Successful
	HYPERSPECTRAL_NO_ERR,
	///Header file not found
	HYPERSPECTRAL_HEADER_FILE_NOT_FOUND,
	///Image file not found
	HYPERSPECTRAL_IMAGE_FILE_NOT_FOUND,
	///Could not find the requested property in the header file
	HYPERSPECTRAL_HDR_PROPERTY_NOT_FOUND,
	///Interleave in header file not supported by this software
	HYPERSPECTRAL_INTERLEAVE_UNSUPPORTED,
	///Datatype in hyperspectral image not supported
	HYPERSPECTRAL_DATATYPE_UNSUPPORTED,
	///mmap failed
	HYPERSPECTRAL_MMAP_FAILED,
	///Failed to read from stdin
	HYPERSPECTRAL_STDIN_FAILED,
	///???
	HYPERSPECTRAL_FILE_READING_ERROR
};

/**
 * Returns error message for specified error type.
 *
 * \param error_code Error code
 * \return Error message
 **/
const char *hyperspectral_error_message(hyperspectral_err_t error_code);

/**
 * Interleave, BIL. BSQ and BIL so far not supported.
 **/
enum interleave_t{BIL_INTERLEAVE, BIP_INTERLEAVE, BSQ_INTERLEAVE};

/**
 * Hyperspectral data types as they are stored within the hyperspectral image file.
 **/
enum datatype_t{HYPERSPECTRAL_DATATYPE_FLOAT = 4, HYPERSPECTRAL_DATATYPE_UINT16_T = 12};

/**
 * Container for hyperspectral header file.
 **/
struct hyspex_header {
	///Interleave of image
	interleave_t interleave;
	///Number of pixels in the across-track axis
	int samples;
	///Number of wavelength bands
	int bands;
	///Number of lines in the image (along-track)
	int lines;
	///Offset of the image data from start of the hyperspectral file
	int offset;
	///Wavelengths
	std::vector<float> wlens;
	///Default (RGB) bands
	std::vector<int> default_bands;
	///Datatype of values in hyperspectral file
	int datatype;
};

/**
 * Get element size represented by the data type specified in the hyperspectral header.
 *
 * \param header Header
 * \return Element size
 **/
size_t hyperspectral_element_size(struct hyspex_header header);

/**
 * Read header information from file. Will read from stdin if filename is empty.
 *
 * \param filename Filename
 * \param header Output header container
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_header(const char *filename, struct hyspex_header *header);

/**
 * Read header information from stdin (for piping hyperspectral data from a process that is using
 * hyperspectral_file_from_stdout() to generate a hyperspectral_file_t instance and writing
 * to that).
 *
 * \param header Output header
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_header_from_stdin(struct hyspex_header *header);

/**
 * M-mapping of hyperspectral file.
 **/
typedef struct {
	///File descriptor
	int file_descriptor;
	///File size that was mapped
	size_t mapped_size;
	///Pointer to mapped memory of file
	void *mmapped_file;
	///Header corresponding to the file
	struct hyspex_header header;
} hyperspectral_image_t;

/**
 * NOTE: The direct access functions below are potentially extremely
 * unsafe, but are part of the public API for performance reasons in
 * specific applications. The safer (but more memory-extensive) way
 * is to use the reading functions that read hyperspectral data
 * directly into float* arrays.
 *
 * No protection or sane access API supplied, FIXME?
 **/

/**
 * Get pointer to specific line in hyperspectral image, given that the data
 * type is float. NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \param line Line
 * \param Pointer to mapped region of hyperspectral image corresponding to the
 * start of the specified line, assuming data type is float.
 **/
float *hyperspectral_float_line(hyperspectral_image_t* data, int line);

/**
 * Get pointer to entire hyperspectral image, given that the data type is
 * float. NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \return Pointer to mmapped region of hyperspectral image corresponding to
 * the start of the image, assuming datatype is float.
 **/
float *hyperspectral_float(hyperspectral_image_t *data);

/**
 * Get pointer to specific line in hyperspectral image, given that the data
 * type is uint16_t. NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \param line Line
 * \param Pointer to mapped region of hyperspectral image corresponding to the
 * start of the specified line, assuming data type is uint16_t.
 **/
uint16_t *hyperspectral_uint16_line(hyperspectral_image_t *data, int line);

/**
 * Get pointer to entire hyperspectral image, given that the data type is
 * uint16_t. NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \return Pointer to mmapped region of hyperspectral image corresponding to
 * the start of the image, assuming datatype is uint16_t.
 **/
uint16_t *hyperspectral_uint16(hyperspectral_image_t *data);

/**
 * Get pointer to entire hyperspectral image as a char array (conversions
 * will have to be done later).
 * NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \return Pointer to mmapped region of hyperspectral image corresponding to
 * the start of the image.
 **/
char *hyperspectral_char(hyperspectral_image_t *data);

/**
 * Get pointer to specific line in hyperspectral image as a char array (conversions
 * will have to be done later).
 * NB: No checks are performed, application will segfault if the
 * returned array is modified, as it is a direct pointer to the write-protected
 * mmapped memory.
 *
 * \param data Hyperspectral data
 * \return Pointer to mmapped region of hyperspectral image corresponding to
 * the start of the image.
 **/
char *hyperspectral_char_line(hyperspectral_image_t *data, int line);

/**
 * mmap hyperspectral file into a hyperspectral_image_t container. Use access functions above.
 *
 * \param filename Input filename
 * \param header Header information
 * \param data Output data
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header header, hyperspectral_image_t *data);

/**
 * mmap hyperspectral file into a hyperspectral_image_t container, and read header in one go.
 *
 * \param filename Input filename
 * \param header Output header information
 * \param data Output data
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_header_and_image(const char *filename, struct hyspex_header *header, hyperspectral_image_t *data);

/**
 * Free mmapped data associated with a hyperspectral_image_t container.
 *
 * \param data Data to free
 **/
void hyperspectral_free_data(hyperspectral_image_t *data);

/**
 * Image subset convenience struct for specifying a subset of the hyperspectral image file for reading.
 **/
struct image_subset {
	///Start pixel in the sample direction
	int start_sample;
	///End pixel in the sample direction
	int end_sample;
	///Start pixel in the line direction
	int start_line;
	///End pixel in the line direction
	int end_line;
	///Start wavelength band
	int start_band;
	///End wavelength band
	int end_band;
};

/**
 * Get image_subset corresponding to the entire image represented by hyspex_header.
 *
 * \param header Hyperspectral header
 * \return Image subset
 **/
struct image_subset hyperspectral_generate_subset(struct hyspex_header header);

/**
 * Get image subset corresponding to the range subset.
 *
 * \param header Hyperspectral header
 * \param start_line Start line
 * \param num_lines Number of lines in subset
 * \param start_sample Start sample
 * \param num_samples Number of samples in subset
 * \param start_band Start band
 * \param num_bands Number of bands in subset
 * \return Image subset
 **/
struct image_subset hyperspectral_generate_subset(struct hyspex_header header, int start_line, int num_lines, int start_sample, int num_samples, int start_band, int num_bands);

/**
 * Construct new header based on subset and original header.
 *
 * \param header Original header
 * \param image_subset Subset
 * \return Subsetted image header
 **/
struct hyspex_header hyperspectral_header_from_subset(struct hyspex_header header, struct image_subset);

/**
 * Get sizes corresponding to the specified image subset.
 *
 * \param image_subset Image subset
 * \param num_lines Returned number of lines
 * \param num_bands Returned number of bands
 * \param num_samples Returned number of samples
 **/
void hyperspectral_get_size(struct image_subset image_subset, int *num_lines, int *num_bands, int *num_samples);

/**
 * Allocate hyperspectral image array as according to sizes in header.
 *
 * \param header Hyperspectral header
 * \return Allocated hyperspectral float array
 **/
float *hyperspectral_alloc_float(struct hyspex_header header);

/**
 * Allocate hyperspectral image array as according to the image subset.
 *
 * \param image_subset Image subset
 * \return Allocated hyperspectral float array
 **/
float *hyperspectral_alloc_float(struct image_subset image_subset);

/**
 * Read hyperspectral image from file.
 *
 * Since image is converted and read into a float array, some memory overhead is to be expected.
 * If memory is an issue, use the API associated with hyperspectral_image_t.
 *
 * \param filename Filename
 * \param header Header, already read from file using hyperspectral_read_header
 * \param image_subset Specified image subset
 * \param data Output data, preallocated to neccessary size (use hyperspectral_image_alloc).
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header *header, struct image_subset subset, float *data);

/**
 * Overloaded version of hyperspectral_read_image where it isn't neccessary to supply subset information, the full image will be read.
 *
 * \param filename Filename
 * \param header Header information
 * \param data Output image data
 * \return HYPERSPECTRAL_NO_ERR on success
 **/
hyperspectral_err_t hyperspectral_read_image(const char *filename, struct hyspex_header *header, float *data);

#include <fstream>

/**
 * Hyperspectral file representation, for writing data to file over time and repeatedly.
 *
 * Preferred way of writing hyperspectral files.
 **/
typedef struct {
	///Output filename
	std::string filename;
	///.img file
	FILE *image_file;
	///Image header
	struct hyspex_header header;
	///Number of written bytes in total
	size_t written_bytes;
} hyperspectral_file_t;

/**
 * Open hyperspectral file for writing. Will overwrite any existing file. Assumes float image.
 *
 * \param filename Filename
 * \param bands Number of bands
 * \param samples Number of samples
 * \param wlens Wavelengths
 * \return File container
 **/
hyperspectral_file_t hyperspectral_open_write_file(std::string filename, int bands, int samples, std::vector<float> wlens);

/**
 * Prepare hyperspectral_file_t for writing to stdout. Will write header information to stdout,
 * and subsequent calls to hyperspectral_write_to_file on this hyperspectral_file_t instance will
 * write to stdout. It is safe to call hyperspectral_close_write_file on this instance. Assumes float image.
 *
 * \param lines Number of lines. Will be used in the header writing, but doesn't really matter unless it matters on the pipe end
 * \param bands Number of bands
 * \param samples Number of samples
 * \param wlens Wavelengths
 * \return File container
 **/
hyperspectral_file_t hyperspectral_open_stdout_file(int lines, int bands, int samples, std::vector<float> wlens);

/**
 * Calls open_write_file if filename is non-empty, otherwise opens stdout. Assumes image to be float.
 *
 * \param filename Filename. Will write to stdout if empty
 * \param lines Expected number of lines
 * \param bands Number of bands
 * \param samples Number of samples
 * \param wlens Wavelengths
 * \return File container
 **/
hyperspectral_file_t hyperspectral_open(std::string filename, int lines, int bands, int samples, std::vector<float> wlens);

/**
 * Open file writer.
 *
 * \param filename Filename. Writes to stdout if empty
 * \param header Stipulated hyperspectral header
 **/
hyperspectral_file_t hyperspectral_open(std::string filename, struct hyspex_header header);

/**
 * Write data to hyperspectral file.
 *
 * \param file File
 * \param num_bytes Number of bytes
 * \param data Binary data
 * \return Number of written bytes
 **/
size_t hyperspectral_write_to_file(hyperspectral_file_t *file, size_t num_bytes, const char *data);

/**
 * Close hyperspectral file and write header file to disk containing the file information.
 *
 * \param file File
 **/
void hyperspectral_close_write_file(hyperspectral_file_t *file);

/**
 * Write header information directly to file.
 *
 * Not preferred, use hyperspectral_file_t API.
 *
 * \param filename Filename
 * \param bands Number of bands
 * \param samples Number of samples (across-track)
 * \param lines Number of lines (along-track)
 * \param wlens Wavelength array
 **/
void hyperspectral_write_header(const char *filename, int bands, int samples, int lines, std::vector<float> wlens);

/**
 * Write hyperspectral float image directly to file, and close the file.
 *
 * Not preferred, use hyperspectral_file_t API.
 *
 * \param filename Filename
 * \param bands Bands
 * \param samples Samples
 * \param lines Lines
 * \param data Image data
 **/
void hyperspectral_write_image(const char *filename, int bands, int samples, int lines, float *data);

void hyperspectral_write_image(std::string filename, struct hyspex_header header, char *data);

void hyperspectral_write_header(std::string filename, struct hyspex_header header);

#endif
