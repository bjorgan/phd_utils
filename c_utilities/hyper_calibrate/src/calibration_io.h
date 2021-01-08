#ifndef CALIBRATION_IO_H_DEFINED
#define CALIBRATION_IO_H_DEFINED

#include "calibration.h"

///Postfix/file extension of filenames for calibration info files
const std::string CALIBRATION_INFO_FILE_POSTFIX = ".calibration_info";

/**
 * Output calibration info to file, or stdout if basefilename is empty.
 *
 * \param calibration Calibration info
 * \param basefilename Output filename
 **/
void hyperspectral_calibration_info_to_file(const hyperspectral_calibration_info_t *calibration, std::string basefilename);


/**
 * Read calibration info from file.
 *
 * \param filename Filename from which to read calibration info.
 * \return Calibration info
 **/
hyperspectral_calibration_info_t hyperspectral_calibration_info_from_file(std::string filename);

/**
 * Check whether calibration info is valid and can be applied on the image represented by the hyperspectral header.
 *
 * \param header Hyperspectral header
 * \param calibration_info Calibration info to check
 * \return True if valid, false otherwise
 **/
bool hyperspectral_calibration_info_valid(struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info);

/**
 * Create calibration info struct.
 *
 * \param num_samples Number of samples
 * \param num_bands Number of bands
 * \param calibration_array Calibration array, is copied to the appropriate field in the struct
 * \param image_subset Image subset from which the calibration info was integrated, if relevant. Enter empty struct if not relevant
 * \return Calibration info
 **/
hyperspectral_calibration_info_t calibration_create_info(int num_samples, int num_bands, double *calibration_array, struct image_subset image_subset);


#endif
