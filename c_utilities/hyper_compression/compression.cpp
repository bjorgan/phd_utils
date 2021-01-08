#include "compression.h"
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <hyperspectral/readimage_stdin.h>

int bil_index(const struct hyspex_header &header, int line, int sample, int band)
{
	return line*header.samples*header.bands + band*header.samples + sample;
}

int bip_index(const struct hyspex_header &header, int line, int sample, int band)
{
	return line*header.samples*header.bands + sample*header.bands + band;
}

void hyperspectral_band(int band, const hyperspectral_compression_t *image, float *ret_image)
{
	float *band_inverse_transform = image->transform->inverse_transform_with_band_selection + band*image->image_header.bands;

	int lines = image->image_header.lines;
	int samples = image->image_header.samples;
	int bands = image->num_bands_in_inverse;

	for (int i=0; i < samples*lines; i++) {
		ret_image[i] = image->transform->means[band];
	}

	cblas_sgemv(CblasRowMajor, CblasNoTrans, samples*lines, bands, 1.0f, image->transformed_image_bip, bands, band_inverse_transform, 1, 1.0f, ret_image, 1);
}

void hyperspectral_pixel(int line, int sample, const hyperspectral_compression_t *image, float *ret_spectrum)
{
	float *z = image->transformed_image_bip + image->num_bands_in_inverse*(line*image->image_header.samples + sample);
	memcpy(ret_spectrum, image->transform->means, image->image_header.bands*sizeof(float));
	cblas_sgemv(CblasRowMajor, CblasNoTrans, image->image_header.bands, image->num_bands_in_inverse, 1.0f, image->transform->inverse_transform_with_band_selection, image->image_header.bands, z, 1, 1.0f, ret_spectrum, 1); //multiply the transform with z, and add the mean (mean contained in ret_spectrum)
}

void mnf_add_mean_to_line(const float *ones_samples, const float *means, int bands, int samples, float *line);

void hyperspectral_line(int line, const hyperspectral_compression_t *image, float *ret_line)
{
	size_t samples_in_line = image->image_header.samples*image->num_bands_in_inverse;

	int bands = image->num_bands_in_inverse;
	int samples = image->image_header.samples;

	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, image->image_header.bands, samples, bands, 1.0f, image->transform->inverse_transform_with_band_selection, image->image_header.bands, image->transformed_image_bip + line*samples_in_line, bands, 0.0f, ret_line, samples);

	mnf_add_mean_to_line(image->transform->ones_samples, image->transform->means, image->image_header.bands, samples, ret_line);
}

#include "progress_bar.h"

hyperspectral_compression_t *allocate_compression(struct hyspex_header header, int num_bands_in_inverse)
{
	hyperspectral_compression_t *compression = new hyperspectral_compression_t;
	compression->num_bands_in_inverse = num_bands_in_inverse;
	compression->transformed_header = header;
	compression->transformed_header.bands = num_bands_in_inverse;
	compression->image_header = header;
	compression->transform = mnf_transform_create(header, num_bands_in_inverse);
	compression->mean_dot_inv = new double[header.bands]();
	compression->mean_euclidean_length_squared = 0;
	compression->transformed_image_bip = new float[header.lines*header.samples*num_bands_in_inverse]();
	return compression;
}

void apply_premultiplication(hyperspectral_compression_t *compression)
{
	struct hyspex_header header = compression->image_header;
	//premultiplication of some vectors for dot product convenience
	for (int i=0; i < header.bands; i++) {
		compression->mean_dot_inv[i] = 0;
		compression->mean_euclidean_length_squared += pow(compression->transform->means[i], 2);
		for (int j=0; j < header.bands; j++) {
			compression->mean_dot_inv[i] += compression->transform->means[j]*compression->transform->inverse_transform[i*header.bands + j];
		}
	}
}

hyperspectral_compression_t* hyperspectral_create_compression(const mnf_statistics_t *statistics, hyperspectral_reader_t *reader, int num_bands_in_inverse, int verbosity)
{
	struct hyspex_header header = reader->header;
	hyperspectral_compression_t *compression = allocate_compression(header, num_bands_in_inverse);
	float *line = new float[header.samples*header.bands]();

	struct hyspex_header line_header = header;
	line_header.lines = 1;

	//prepare MNF transform
	mnf_transform_calculate(header, statistics, compression->transform);
	mnf_orthonormalize_transform(compression->num_bands_in_inverse, header, compression->transform);
	apply_premultiplication(compression);

	//transform image, but keep only first relevant bands
	for (int i=0; i < header.lines; i++) {
		if (verbosity > 0) print_progress_bar("Compress image.", i-1, i, header.lines-1);
		hyperspectral_line_data_float(reader, i, line);
		mnf_run_forward(compression->transform, line_header, line);

		//copy to BIP array
		for (int j=0; j < header.samples; j++) {
			for (int k=0; k < num_bands_in_inverse; k++) {
				compression->transformed_image_bip[i*num_bands_in_inverse*header.samples + j*num_bands_in_inverse + k] = line[k*header.samples + j];
			}
		}
	}

	delete [] line;

	return compression;
}

hyperspectral_compression_t* hyperspectral_create_compression(hyperspectral_reader_t *reader, int num_bands_in_inverse, int verbosity)
{
	//prepare MNF statistics
	struct hyspex_header header = reader->header;
	size_t element_size = hyperspectral_element_size(header);
	size_t line_size = header.samples*header.bands*element_size;
	mnf_statistics_t *statistics = mnf_statistics_create(header);
	float *line = new float[header.samples*header.bands]();
	struct hyspex_header line_header = header;
	line_header.lines = 1;
	for (int i=0; i < header.lines; i++) {
		if (verbosity > 0) print_progress_bar("Calculate MNF transform.", i-1, i, header.lines-1);
		hyperspectral_line_data_float(reader, i, line);
		mnf_update_statistics(statistics, line_header, line);
	}
	delete [] line;

	hyperspectral_compression_t *ret_compression = hyperspectral_create_compression(statistics, reader, num_bands_in_inverse, verbosity);
	mnf_free_statistics(&statistics);
	return ret_compression;
}

hyperspectral_compression_t* hyperspectral_create_compression(struct hyspex_header header, const char *image_data, int num_bands_in_inverse, int verbosity)
{
	//create hyperspectral reader where we pretend that we have already read image_data from stdin.
	hyperspectral_reader_t reader;
	hyperspectral_stdin_t dummy_stdin;
	dummy_stdin.curr_line_in_stdin = header.lines;
	dummy_stdin.stream_closed = true;
	dummy_stdin.should_accumulate_data = true;
	dummy_stdin.header = header;
	dummy_stdin.accumulated_data = (char*)image_data;
	reader.stdin_stream = &dummy_stdin;
	reader.read_from_stdin = true;
	reader.header = header;

	return hyperspectral_create_compression(&reader, num_bands_in_inverse, verbosity);
}

void hyperspectral_free_compression(hyperspectral_compression_t **compressed_image)
{
	mnf_free_transform(&((*compressed_image)->transform));
	delete [] (*compressed_image)->transformed_image_bip;
	delete [] (*compressed_image)->mean_dot_inv;
}

double hyperspectral_dot_product(const hyperspectral_compression_t *image, coordinate i, coordinate j)
{
	double retval = image->mean_euclidean_length_squared;
	for (int k=0; k < image->num_bands_in_inverse; k++) {
		size_t ind_i = bip_index(image->transformed_header, i.first, i.second, k);
		size_t ind_j = bip_index(image->transformed_header, j.first, j.second, k);
		float z_i = image->transformed_image_bip[ind_i];
		float z_j = image->transformed_image_bip[ind_j];
		retval += z_i*z_j;
		retval += image->mean_dot_inv[k]*(z_i + z_j);
	}
	return retval;
}

double hyperspectral_euclidean_distance_squared(const hyperspectral_compression_t *image, coordinate i, coordinate j)
{
	double retval = 0.0;
	for (int k=0; k < image->num_bands_in_inverse; k++) {
		size_t ind_i = bip_index(image->transformed_header, i.first, i.second, k);
		size_t ind_j = bip_index(image->transformed_header, j.first, j.second, k);
		retval += pow(image->transformed_image_bip[ind_i] - image->transformed_image_bip[ind_j], 2);
	}
	return retval;

}

#include <limits>

void mnf_transform_to_file(std::string filename, struct hyspex_header header, mnf_transform_t *transform)
{
	std::ofstream out_file(filename.c_str());
	out_file.precision(std::numeric_limits<double>::digits10 + 1);

	out_file << "MNF transform";
	out_file << std::endl; for (int i=0; i < header.bands*header.bands; i++) out_file << transform->forward_transform[i] << " ";
	out_file << std::endl; for (int i=0; i < header.bands*header.bands; i++) out_file << transform->inverse_transform[i] << " ";
	out_file << std::endl; for (int i=0; i < header.bands*header.bands; i++) out_file << transform->inverse_transform_with_band_selection[i] << " ";
	out_file << std::endl; for (int i=0; i < header.bands*header.bands; i++) out_file << transform->denoising_transform[i] << " ";
	out_file << std::endl; for (int i=0; i < header.bands; i++) out_file << transform->eigenvalues[i] << " ";
	out_file << std::endl; for (int i=0; i < header.bands; i++) out_file << transform->means[i] << " ";
}

void read_float_line(std::ifstream *file, int samples, float *output_array)
{
	std::string line;
	getline(*file, line);
	std::stringstream ss(line);
	for (int i=0; i < samples; i++) {
		ss >> output_array[i];
	}
}

void mnf_transform_from_file(std::string filename, struct hyspex_header header, mnf_transform_t *transform)
{
	std::ifstream in_file(filename.c_str());
	std::string header_line;
	std::getline(in_file, header_line); //first line is text "MNF transform"
	read_float_line(&in_file, header.bands*header.bands, transform->forward_transform);
	read_float_line(&in_file, header.bands*header.bands, transform->inverse_transform);
	read_float_line(&in_file, header.bands*header.bands, transform->inverse_transform_with_band_selection);
	read_float_line(&in_file, header.bands*header.bands, transform->denoising_transform);
	read_float_line(&in_file, header.bands, transform->eigenvalues);
	read_float_line(&in_file, header.bands, transform->means);

}

#include <unistd.h>
#include <string.h>
#include <libtar.h>
#include <fcntl.h>

#include <dirent.h>

/**
 * Remove temporary (actually any) directory and its contents.
 *
 * \param tmp_dir Path to temporary directory
 **/
void rm_temp_dir(std::string tmp_dir)
{
	DIR *d;
	struct dirent *dir;
	d = opendir(tmp_dir.c_str());
	std::vector<std::string> filenames;
	if (d) {
		while ((dir = readdir(d)) != NULL) {
			if (dir->d_type == DT_REG) {
				filenames.push_back(tmp_dir + "/" + std::string(dir->d_name));
			}
		}
		closedir(d);
	}

	for (int i=0; i < filenames.size(); i++) {
		unlink(filenames[i].c_str());
	}
	rmdir(tmp_dir.c_str());
}

const std::string COMPRESSED_IMAGE_HEADER_BASE = "image_header";
const std::string COMPRESSED_IMAGE_HEADER_FNAME = COMPRESSED_IMAGE_HEADER_BASE + HEADER_EXTENSION;
const std::string COMPRESSED_TRANSFORMED_IMAGE_BASE = "transformed_image";
const std::string COMPRESSED_TRANSFORMED_IMAGE_FNAME = COMPRESSED_TRANSFORMED_IMAGE_BASE + IMAGE_EXTENSION;
const std::string COMPRESSED_TRANSFORMED_HEADER_FNAME = COMPRESSED_TRANSFORMED_IMAGE_BASE + HEADER_EXTENSION;
const std::string COMPRESSED_TRANSFORM_FNAME = "transform.dat";

void hyperspectral_compression_to_file(std::string filename, const hyperspectral_compression_t *compressed_image)
{
	//write temporary files
	char dir_template[] = "/tmp/writecompr.XXXXXX";
	char *tmp_dir = strdup(mkdtemp(dir_template));

	std::string image_header = std::string(tmp_dir) + "/" + COMPRESSED_IMAGE_HEADER_BASE;

	std::string transformed_image = std::string(tmp_dir) + "/" + COMPRESSED_TRANSFORMED_IMAGE_BASE;

	std::string transform = std::string(tmp_dir) + "/" + COMPRESSED_TRANSFORM_FNAME;

	hyperspectral_write_header(image_header, compressed_image->image_header);
	hyperspectral_write_header(transformed_image.c_str(), compressed_image->transformed_header);
	hyperspectral_write_image(transformed_image.c_str(), compressed_image->transformed_header.bands, compressed_image->transformed_header.samples, compressed_image->transformed_header.lines, compressed_image->transformed_image_bip);
	mnf_transform_to_file(transform, compressed_image->image_header, compressed_image->transform);

	image_header += HEADER_EXTENSION;
	std::string transformed_header = transformed_image + HEADER_EXTENSION;
	transformed_image += IMAGE_EXTENSION;

	//create tar archive out of the files
	TAR *tarfile;
	tar_open(&tarfile, filename.c_str(), NULL, O_WRONLY | O_CREAT, 0644, TAR_GNU);
	tar_append_file(tarfile, image_header.c_str(), COMPRESSED_IMAGE_HEADER_FNAME.c_str());
	tar_append_file(tarfile, transformed_header.c_str(), COMPRESSED_TRANSFORMED_HEADER_FNAME.c_str());
	tar_append_file(tarfile, transformed_image.c_str(), COMPRESSED_TRANSFORMED_IMAGE_FNAME.c_str());
	tar_append_file(tarfile, transform.c_str(), COMPRESSED_TRANSFORM_FNAME.c_str());
	tar_append_eof(tarfile);
	tar_close(tarfile);

	//clean up temporary files
	rm_temp_dir(tmp_dir);
	free(tmp_dir);
}

hyperspectral_err_t hyperspectral_read_image_sequential(const char *filename, struct hyspex_header *header, float *data);

hyperspectral_compression_t* hyperspectral_compression_from_file(std::string filename)
{
	//extract tar archive into temporary directory
	TAR *tar;
	if (tar_open(&tar, filename.c_str(), NULL, O_RDONLY, 0, 0) != 0) {
		return NULL;
	}

	char dir_template[] = "/tmp/readcompr.XXXXXX";
	char *tmp_dir_c_str = strdup(mkdtemp(dir_template));
	std::string tmp_dir = std::string(tmp_dir_c_str) + "/";

	if (tar_extract_all(tar, tmp_dir_c_str) != 0) {
		return NULL;
	}

	tar_close(tar);
	free(tmp_dir_c_str);

	//read original image header
	struct hyspex_header header;
	hyperspectral_err_t errcode = hyperspectral_read_header((tmp_dir + COMPRESSED_IMAGE_HEADER_FNAME).c_str(), &header);

	//read number of bands in inverse
	struct hyspex_header transformed_header;
	errcode = hyperspectral_read_header((tmp_dir + COMPRESSED_TRANSFORMED_HEADER_FNAME).c_str(), &transformed_header);
	int num_bands_in_inverse = transformed_header.bands;

	//prepare compression and read image data
	hyperspectral_compression_t *compressed_image = allocate_compression(header, num_bands_in_inverse);
	errcode = hyperspectral_read_image_sequential((tmp_dir + COMPRESSED_TRANSFORMED_IMAGE_FNAME).c_str(), &transformed_header, compressed_image->transformed_image_bip);

	compressed_image->image_header = header;
	compressed_image->transformed_header = transformed_header;

	//read transform
	mnf_transform_from_file(tmp_dir + COMPRESSED_TRANSFORM_FNAME, header, compressed_image->transform);
	apply_premultiplication(compressed_image);

	//remove temporary directory
	rm_temp_dir(tmp_dir);
	return compressed_image;
}

#include <magic.h>

const std::string TAR_MIME_TYPE = "application/x-tar";

bool hyperspectral_file_is_compressed(std::string filepath)
{
	magic_t magic = magic_open(MAGIC_MIME_TYPE);
	magic_load(magic, NULL);
	magic_compile(magic, NULL);
	std::string mime_type = std::string(magic_file(magic, filepath.c_str()));
	if (mime_type == TAR_MIME_TYPE) {
		return true;
	}
	return false;
}
