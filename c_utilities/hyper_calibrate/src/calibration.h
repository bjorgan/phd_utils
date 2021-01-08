#ifndef CALIBRATION_H_DEFINED
#define CALIBRATION_H_DEFINED

#include <hyperspectral/readimage.h>

/**
 * Information about integrated reflectance standard.
 **/
typedef struct {
 	///Integrated calibration values over a reflectance standard. Double is used, but the difference between float and double is minuscle.
	double *calibration_array;
	///Number of samples present in calibration_array. Should either correspond to width of image or be exactly 1.
	int num_samples;
	///Number of bands in calibration_array. Should correspond to the number of bands in the image.
	int num_bands;
	///Whether calibration_array has been adjusted against true reflectance standard specs.
	bool adjusted;
	///Image subset used for reflectance standard, if relevant
	struct image_subset refl_standard_positions;
} hyperspectral_calibration_info_t;

/**
 * Free memory associated with calibration_info.
 *
 * \param calibration_info Calibration info
 **/
void hyperspectral_calibration_info_free(hyperspectral_calibration_info_t *calibration_info);

/**
 * Used for specifying direction for which we try to search for the reflectance standard.
 **/
enum search_direction{
	START_FROM_END,
	START_FROM_START
};

/**
 * Find a calibration standard that lies in either extreme end of the image, and is covering the entire field of view.
 *
 * Will not try to search for the start of the reflectance standard, only the end. If there are artifacts
 * in the image or changes at the end or start before or after the reflectance standard, this code
 * will not work properly.
 *
 * \param direction Direction from which to start looking (end or start of image)
 * \param header Image header
 * \param image_data Image data (BIL-interleaved)
 * \param verbosity Verbosity level
 * \return Positions for reflectance standard
 **/
struct image_subset calibration_find_FOV_refl_standard(enum search_direction direction, struct hyspex_header header, const char *image_data, int verbosity = 0);

/**
 * Modify reflectance standard ROI in sample direction, by applying the same algorithm for finding start and end lines in across-track direction.
 *
 * \param image_subset Image subset found using calibration_find_FOV_refl_standard()
 * \param header Image header
 * \param image_data Image data (BIL-interleaved)
 * \param verbosity Verbosity level
 **/
void calibration_crop_refl_standard_in_sample_direction(struct image_subset &image_subset, struct hyspex_header header, const char *image_data, int verbosity);

/**
 * Integrate an input image in order to find the corresponding calibration array.
 * \param header Image header
 * \param image_data Image data
 * \param subset Image subset we will integrate over
 * \param verbosity Verbosity level
 * \return Calibration info
 **/
hyperspectral_calibration_info_t calibration_integrate(struct hyspex_header header, const char *image_data, struct image_subset subset, int verbosity = 0);

/**
 * Integrate an input image in order to find the corresponding calibration array. Integrate over full subset, and use same calibration spectrum in all pixels of output array.
 * \param header Image header
 * \param image_data Image data
 * \param subset Image subset we will integrate over
 * \param verbosity Verbosity level
 * \return Calibration info
 **/
hyperspectral_calibration_info_t calibration_integrate_all(struct hyspex_header header, const char *image_data, struct image_subset subset, int verbosity);

/**
 * Adjust calibration array against specification file for the reflectance standard.
 * \param specs_filename Path to reflectance standard specifications file
 * \param header Image header
 * \param calibration_info Calibration info to be adjusted
 **/
void calibration_adjust(std::string specs_filename, struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info);

/**
 * Adjust calibration array against a constant specification value for the reflectance standard.
 * \param specification_value Value corresponding to the specification of the standard
 * \param header Image header
 * \param calibration_info Calibration info to be adjusted
 **/
void calibration_adjust_with_constant_value(double specification_value, struct hyspex_header header, hyperspectral_calibration_info_t *calibration_info);


/**
 * Calibrate single line of data according to input calibration array.
 *
 * \param calibration_info Calibration info
 * \param header Image header
 * \param input_data Input data
 * \param output_data Output calibrated data
 **/
void calibration_run_single_line(const hyperspectral_calibration_info_t *calibration_info, struct hyspex_header header, const char *input_data, float *output_data);

#ifdef WITH_OPENCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

/**
 * Generate rgb image where borders of reflectance standard are indicated, for validication.
 * \param header Image header
 * \param rgb_data BSQ-interleaved RGB data
 * \param subset Reflectance standard position
 * \return RGB image where borders of the reflectance standard are marked with a red line.
 **/
cv::Mat calibration_generate_rgb_image(struct hyspex_header header, float *rgb_data, struct image_subset subset);
#endif

#endif
