#ifndef UV_CAMERA_H_DEFINED
#define UV_CAMERA_H_DEFINED

#include <Jai_Factory.h>
#include <vector>

/**
 * Get camera frame size.
 **/
SIZE uv_frame_size(CAM_HANDLE *camera);

/**
 * Close factory and camera handles.
 *
 * \param factory Factory handle, to be closed
 * \param camera Camera handle, to be closed
 **/
void uv_close_all(FACTORY_HANDLE *factory, CAM_HANDLE *camera);

/**
 * Get min/max ranges for absolute exposure times in us.
 *
 * \param camera Camera handle
 * \param min Returned minimum in us
 * \param max Returned maximum in us
 **/
void uv_camera_exposure_time_minmax(CAM_HANDLE *camera, double *min, double *max);

/**
 * Get min/max ranges for raw exposure times.
 *
 * \param camera Camera handle
 * \param min Returned raw minimum
 * \param max Returned raw maximum
 **/
void uv_camera_raw_exposure_time_minmax(CAM_HANDLE *camera, long *min, long *max);

/**
 * Get exposure time in absolute units.
 *
 * \param camera Camera handle
 * \return Exposure time in microseconds
 **/
double uv_camera_exposure_time(CAM_HANDLE *camera);

/**
 * Get raw exposure time.
 *
 * \param camera Camera handle
 * \return Raw exposure time
 **/
long uv_camera_raw_exposure_time(CAM_HANDLE *camera);

/**
 * Set exposure time in absolute units.
 *
 * \param camera Camera handle
 * \param exposure_time Exposure time in microseconds.
 **/
void uv_camera_set_exposure_time(CAM_HANDLE *camera, double exposure_time);

/**
 * Get min/max ranges for possible gains.
 *
 * \param camera Camera handle
 * \param min Returned minimum
 * \param max Returned maximum
 **/
void uv_camera_gain_minmax(CAM_HANDLE *camera, long *min, long *max);

/**
 * Set gain in absolute units.
 *
 * \param camera Camera handle
 * \param input_gain Gain
 **/
void uv_camera_set_gain(CAM_HANDLE *camera, long input_gain);

/**
 * Get gain in absolute units.
 *
 * \param camera Camera handle
 * \return Gain
 **/
long uv_camera_gain(CAM_HANDLE *camera);

/**
 * Send start acquisition command to camera.
 *
 * \param camera Camera handle
 **/
void uv_camera_start_acquisition(CAM_HANDLE *camera);

/**
 * Send stop acquisition command to camera.
 *
 * \param camera Camera handle
 **/
void uv_camera_stop_acquisition(CAM_HANDLE *camera);

///Length of histogram array
const int HISTOGRAM_LENGTH = 256;

/**
 * Histogram structure.
 **/
typedef struct {
	///Pixel count for each value from 0 to 255
	int counts[HISTOGRAM_LENGTH];
} histogram_mono8_t;

/**
 * Calculate histogram from image. Assumes monochromatic 8 bit image.
 *
 * \param image Converted image
 * \param histogram Returned histogram
 **/
void uv_camera_histogram(J_tIMAGE_INFO *image, histogram_mono8_t *histogram);

/**
 * Copy image data from src to dst.
 **/
void copy_image(J_tIMAGE_INFO *dst, J_tIMAGE_INFO *src);

/**
 * Setting for verbosity.
 **/
enum verbosity_t {
	///No output
	NO_TERMINAL_OUTPUT = 0,
	///Verbose output
	VERBOSE_OUTPUT = 1
};

/**
 * Save list of images to file. If list contains only one element, the supplied filename is used, while if list contains
 * multiple elements, the filename is postfixed with a frame number.
 * The filename is always postfixed with ".tiff".
 *
 * \param images List of images
 * \param filename Filename for saving
 * \param verbose Whether to print verbose output to terminal
 **/
void uv_save_images(std::vector<J_tIMAGE_INFO> images, std::string filename, verbosity_t verbosity=NO_TERMINAL_OUTPUT);

/**
 * Postfix a string with current timestamp and various camera parameters.
 *
 * \param base_string Base string
 * \param camera Camera handle
 * \return String postfixed with various information
 **/
std::string uv_camera_add_timestamps_and_info(CAM_HANDLE *camera, std::string base_string);

#endif
