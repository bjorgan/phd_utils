#include <hyperspectral/mnf.h>
#include <hyperspectral/readimage.h>
#include <hyperspectral/readimage_stdin.h>
#include <utility>

/**
 * Compressed hyperspectral image
 **/
typedef struct {
	///Number of bands kept in the transformed image
	int num_bands_in_inverse;
	///Transformed, compressed image, containing only `num_bands_in_inverse` bands. BIP-interleaved, for convenience,
	//though all backtransformations will be BIL.
	float *transformed_image_bip;
	///Header of original image
	struct hyspex_header image_header;
	///Header of transformed image, copied from original image but with number of bands set to `num_bands_in_inverse`
	struct hyspex_header transformed_header;
	///MNF transform used to compress the image. The first `num_bands_in_inverse` bands in the inverse transform are orthonormalized, and the forward transform adjusted accordingly.
	mnf_transform_t *transform;
	///Dot product between the mean and each column vector in the inverse transform
	double *mean_dot_inv;
	///Euclidean length of the mean vector
	double mean_euclidean_length_squared;
} hyperspectral_compression_t;

/**
 * Unpack band image.
 *
 * \param band Band index
 * \param image Compressed image
 * \param ret_image Returned image array (header.lines X header.samples)
 **/
void hyperspectral_band(int band, const hyperspectral_compression_t *image, float *ret_image);

/**
 * Unpack pixel.
 *
 * \param line Line coordinate
 * \param sample Sample coordinate
 * \param image Compressed image
 * \param ret_spectrum Returned pixel spectrum (header.bands)
 **/
void hyperspectral_pixel(int line, int sample, const hyperspectral_compression_t *image, float *ret_spectrum);

/**
 * Unpack line.
 *
 * \param line Line coordinate
 * \param image Compressed image
 * \param ret_line Returned hyperspectral line (header.bands X header.samples)
 **/
void hyperspectral_line(int line, const hyperspectral_compression_t *image, float *ret_line);

/**
 * Create compression out of input hyperspectral image.
 *
 * \param header Hyperspectral image header
 * \param image_data Hyperspectral image data
 * \param num_bands_in_inverse Number of bands to use of the MNF transformed image
 * \return Compressed image
 **/
hyperspectral_compression_t* hyperspectral_create_compression(struct hyspex_header header, const char *image_data, int num_bands_in_inverse, int verbosity = 0);

hyperspectral_compression_t* hyperspectral_create_compression(const mnf_statistics_t *statistics, hyperspectral_reader_t *reader, int num_bands_in_inverse, int verbosity = 0);

hyperspectral_compression_t* hyperspectral_create_compression(hyperspectral_reader_t *reader, int num_bands_in_inverse, int verbosity = 0);

/**
 * Free memory associated with compressed image.
 *
 * \param compressed_image Compressed image to free
 **/
void hyperspectral_free_compression(hyperspectral_compression_t **compressed_image);

/**
 * Write hyperspectral compression to file.
 *
 * \param filename Output filename
 * \param compressed_image Compressed image
 **/
void hyperspectral_compression_to_file(std::string filename, const hyperspectral_compression_t *compressed_image);
hyperspectral_compression_t* hyperspectral_compression_from_file(std::string filename);

typedef std::pair<int, int> coordinate;

/**
 * Calculate dot product between two hyperspectral pixels.
 *
 * \param image Compressed image
 * \param i First coordinate
 * \param j Second coordinate
 **/
double hyperspectral_dot_product(const hyperspectral_compression_t *image, coordinate i, coordinate j);

/**
 * Calculate euclidean distance between two hyperspectral pixels.
 *
 * \param image Compressed image
 * \param i First coordinate
 * \param j Second coordinate
 **/
double hyperspectral_euclidean_distance_squared(const hyperspectral_compression_t *image, coordinate i, coordinate j);

/**
 * Returns true if specified file path points to a compressed file.
 *
 * \param filepath Path to file
 * \return True if file is compressed and readable by this library
 **/
bool hyperspectral_file_is_compressed(std::string filepath);
