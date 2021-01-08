/**
 * Asgeir Bjorgan, NTNU, 2018.
 *
 * For accessing the JAI camera SDK for acquiring images using a Jai CM-140GE-UV camera.
 *
 * The SDK manual specifies the general functions, but does not say anything
 * about the parameter names, as they are defined by the camera itself.  Use
 * the JAI Control tool to see tunable parameters, and what names should be
 * used (actually specified by the GigE-standard?).
 **/

#include "uv_camera.h"
#include <stdio.h>
#include <conio.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include <getopt.h>

#include "camera_control_widget.h"
#include <QApplication>

#include <ctime>

///Prefix in the device strings that correspond to a filter driver interface.
#define FD_DEVICE_PREFIX "TL=>GevTL , INT=>FD"

///Command for starting acquisition
int8_t *CMD_ACQUISITION_START = (int8_t*)"AcquisitionStart";
///Command for stopping acquisition
int8_t *CMD_ACQUISITION_STOP = (int8_t*)"AcquisitionStop";
///Width property
int8_t *PROP_WIDTH = (int8_t*)"Width";
///Height property
int8_t *PROP_HEIGHT = (int8_t*)"Height";
///Trigger mode property
int8_t *PROP_TRIGGER_MODE = (int8_t*)"TriggerMode";
///Off value for trigger mode
int8_t *TRIGGER_MODE_OFF = (int8_t*)"Off";
///Exposure time property node name
int8_t *NODE_EXPOSURE_TIME = (int8_t*)"ExposureTimeAbs";
///Raw exposure time property node name
int8_t *NODE_RAW_EXPOSURE_TIME = (int8_t*)"ExposureTimeRaw";
///Gain property node name
int8_t *NODE_GAIN = (int8_t*)"GainRaw";
///Pixel format property node name
int8_t *NODE_PIXEL_FORMAT = (int8_t*)"PixelFormat";
///Monochrome 8bit pixelformat
int8_t *PIXELFORMAT_MONO_8 = (int8_t*)"Mono8";
///Monochrome 10bit packed pixelformat
int8_t *PIXELFORMAT_MONO_10_PACKED = (int8_t*)"Mono10Packed";

/**
 * Prepare JAI Factory. Will call exit(1) on failure.
 *
 * \param factory Factory, to be initialized
 **/
void uv_prepare_factory(FACTORY_HANDLE *factory);

/**
 * Connect to the first available camera having a filter driver interface. Will call exit(1) on failure.
 *
 * Each camera has a socket driver and a filter driver, but filter driver is according to the specifications
 * the most efficient way of communicating with the camera.
 *
 * \param factory Factory, as initialized using uv_prepare_factory()
 * \param camera Camera to be initialized, functioning as a handle for the connected camera
 **/
void uv_prepare_camera(FACTORY_HANDLE *factory, CAM_HANDLE *camera);

/**
 * Used as a callback object for JAI's stream function. Collects camera frames received through the callback function,
 * for acquisition and saving. The live stream in uv_live_view uses a different mechanism.
 **/
class image_collector {
	public:
		/**
		 * Constructor.
		 **/
		image_collector();

		/**
		 * Callback function, used for receiving frames from the JAI camera stream.
		 *
		 * \param raw_image_info Raw image, will be converted to image and stored in converted_frames
		 **/
		void stream_callback(J_tIMAGE_INFO *raw_image_info);

		/**
		 * Get current number of stored frames.
		 **/
		size_t get_num_frames(){return num_frames;};

		/**
		 * Get all stored frames.
		 **/
		std::vector<J_tIMAGE_INFO> get_images(){return converted_frames;};
	private:
		///Current number of frames
		size_t num_frames;
		///Converted and stored frames
		std::vector<J_tIMAGE_INFO> converted_frames;
};

/**
 * Get node corresponding to exposure time settings. Can issue commands to this node in order to get and set settings.
 *
 * \param camera Camera handle
 * \return Node handle for exposure time in absolute units (microseconds)
 **/
NODE_HANDLE uv_camera_exposure_time_property_node(CAM_HANDLE *camera);

/**
 * Get node corresponding to gain settings. Can issue commands to this node in order to get and set settings.
 *
 * \param camera Camera handle
 * \return Node handle for gain
 **/
NODE_HANDLE uv_camera_gain_property_node(CAM_HANDLE *camera);

/**
 * Set pixel format (e.g. monochrome 8 bit, 10 bit, ...).
 *
 * \param camera Camera handle
 * \param pixel_format_string Pixel format string, see JAI Control tool for which strings are possible (not defined in manual)
 **/
void uv_camera_set_pixel_format(CAM_HANDLE *camera, int8_t *pixel_format_string);

/**
 * Set pixel format of camera to be monochrome 8 bit (0-255).
 *
 * \param camera Camera handle
 **/
void uv_camera_set_mono8bit(CAM_HANDLE *camera);

/**
 * Set pixel format of camera to be monochrome 10 bit packed (will be 16 bit in output files, so 0 to something high).
 *
 * \param camera Camera handle
 **/
void uv_camera_set_mono10bit(CAM_HANDLE *camera);

/**
 * Get current pixel format of camera.
 *
 * \param camera Camera handle
 * \return Pixel format
 **/
std::string uv_camera_pixel_format(CAM_HANDLE *camera);

/**
 * Acquire images from the camera.
 *
 * \param camera Camera handle
 * \param num_frames Number of frames to collect
 * \return Vector over acquired and converted frames from the camera
 **/
std::vector<J_tIMAGE_INFO> uv_acquire_images(CAM_HANDLE *camera, size_t num_frames = 1);

/**
 * Save JAI image to file.
 *
 * \param image Converted image
 * \param filename Filename
 **/
void uv_save_image(J_tIMAGE_INFO *image, std::string filename);

///Default exposure time
const double DEFAULT_EXPOSURE_TIME = 118.0;
///Default gain setting
const long DEFAULT_GAIN = 0;

/**
 * Show help.
 *
 * \param program_name Executable name
 **/
void show_help(const char *program_name);

///Live view option number
#define OPT_LIVE_VIEW 201
///No acquisition option number
#define OPT_NO_ACQUISITION 202
///Number of frames option number
#define OPT_NUM_FRAMES 203
///Exposure time option number
#define OPT_EXPOSURE_TIME 204
///Gain option number
#define OPT_GAIN 205
///Show parameters option number
#define OPT_SHOW_PARAMETERS 206

int main(int argc, char *argv[])
{
	//command line options
	struct option long_options[] = {
		{"live-view",			no_argument,		0,	OPT_LIVE_VIEW},
		{"no-acquisition",		no_argument,		0,	OPT_NO_ACQUISITION},
		{"show-parameters",		no_argument,		0,	OPT_SHOW_PARAMETERS},
		{"num-frames",			required_argument,		0,	OPT_NUM_FRAMES},
		{"exposure-time",		required_argument,		0,	OPT_EXPOSURE_TIME},
		{"gain",			required_argument,		0,	OPT_GAIN},
		{"help",			no_argument,		0,	'h'},
		{0, 0, 0, 0}
	};
	char short_options[] = "h";

	//whether to show live view
	bool show_live_view = false;

	//whether to acquire images
	bool acquire_image = true;

	//number of frames to acquire
	long num_frames = 1;

	//exposure time
	double exposure_time = DEFAULT_EXPOSURE_TIME;

	//gain
	double gain = DEFAULT_GAIN;

	//parse command line options
	while (1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, short_options, long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case OPT_LIVE_VIEW: //live view
				show_live_view = true;
				acquire_image = false;
				break;
			case OPT_SHOW_PARAMETERS: //show parameters
			case OPT_NO_ACQUISITION: //no acquisition
				show_live_view = false;
				acquire_image = false;
				break;
			case OPT_NUM_FRAMES: //number of frames
				num_frames = strtod(optarg, NULL);
				break;
			case OPT_EXPOSURE_TIME: //exposure time
				exposure_time = strtod(optarg, NULL);
				break;
			case OPT_GAIN: //gain
				gain = strtod(optarg, NULL);
				break;
			case 'h': //help
				show_help(argv[0]);
				return 0;
				break;
			default:
				exit(1);
		}
	}

	//output filename
	std::string filename;
	if (optind < argc) {
		filename = std::string(argv[optind]);
	}

	if ((filename.size() == 0) && acquire_image) {
		printf("No filename specified: No images will be acquired.\n");
		acquire_image = false;
	}

	//connect to camera (will exit program on failure from within the functions)
	FACTORY_HANDLE factory;
	uv_prepare_factory(&factory);

	CAM_HANDLE camera;
	uv_prepare_camera(&factory, &camera);

	//set parameters
	uv_camera_set_gain(&camera, gain);
	uv_camera_set_exposure_time(&camera, exposure_time);
	uv_camera_set_mono8bit(&camera);

	//show live preview
	if (show_live_view) {
		//display Qt widget for control, will take care of live view starting, stopping and camera destroying
		QApplication qapp(argc, argv);
		camera_control control(&factory, &camera);
		control.show();
		return qapp.exec();
	}

	//acquire images
	if (acquire_image) {
		filename = uv_camera_add_timestamps_and_info(&camera, filename);
		std::vector<J_tIMAGE_INFO> images = uv_acquire_images(&camera, num_frames);

		printf("%ld frames acquired\n", num_frames);
		uv_save_images(images, filename, VERBOSE_OUTPUT);
	}

	//show camera parameters
	if (!acquire_image && !show_live_view) {
		long min_gain, max_gain;
		uv_camera_gain_minmax(&camera, &min_gain, &max_gain);
		printf("Gain range: min = %ld, max = %ld.\n", min_gain, max_gain);
		printf("Current gain: %ld\n", uv_camera_gain(&camera));

		double min_exposure, max_exposure;
		uv_camera_exposure_time_minmax(&camera, &min_exposure, &max_exposure);
		printf("Exposure time range: min = %f us, max = %f us.\n", min_exposure, max_exposure);
		printf("Current exposure time: %f us\n", uv_camera_exposure_time(&camera));

		long min_exposure_raw, max_exposure_raw;
		uv_camera_raw_exposure_time_minmax(&camera, &min_exposure_raw, &max_exposure_raw);
		printf("Raw exposure time range: min = %ld, max = %ld\n", min_exposure_raw, max_exposure_raw);

		printf("Current pixel format: %s\n", uv_camera_pixel_format(&camera).c_str());
	}

	uv_close_all(&factory, &camera);
	return 0;
}

SIZE uv_frame_size(CAM_HANDLE *camera)
{
	SIZE size = {100, 100};
	int64_t int64Val;
	J_Camera_GetValueInt64(*camera, PROP_WIDTH, &int64Val);
	size.cx = (LONG)int64Val;

	J_Camera_GetValueInt64(*camera, PROP_HEIGHT, &int64Val);
	size.cy = (LONG)int64Val;
	return size;
}

void uv_prepare_factory(FACTORY_HANDLE *factory)
{
	J_STATUS_TYPE retval = J_Factory_Open((int8_t*)"" , factory);
	if (retval != J_ST_SUCCESS) {
		fprintf(stderr, "Error accessing JAI SDK factory\n");
		exit(1);
	}
}

void uv_prepare_camera(FACTORY_HANDLE *factory, CAM_HANDLE *camera)
{
	//update camera list
	bool8_t list_has_changed;
	J_STATUS_TYPE retval = J_Factory_UpdateCameraList(*factory, &list_has_changed);
	if (retval != J_ST_SUCCESS) {
		printf("Error accessing camera list\n");
		exit(1);
	}

	//get number of cameras
	uint32_t num_devices;
	retval = J_Factory_GetNumOfCameras(*factory, &num_devices);
	if (retval != J_ST_SUCCESS) {
		printf("Could not get number of cameras.\n");
		exit(1);
	}

	if (num_devices == 0) {
		printf("No cameras found\n");
		exit(1);
	}

	//get the ID of the first available filter driver interface (assume only one camera on the system)
	//(filter driver (FD) is supposed to have better performance than socket driver (SD), so ignoring any SD interfaces)
	int8_t camera_id[J_CAMERA_ID_SIZE] = {0};
	bool found_fd = false;
	for (uint32_t i=0; i < num_devices; i++) {
		uint32_t size = J_CAMERA_ID_SIZE;
		camera_id[0] = 0;
		retval = J_Factory_GetCameraIDByIndex(*factory, i, camera_id, &size);
		if (retval != J_ST_SUCCESS) {
			printf("Failed to get camera ID\n");
			exit(1);
		}

		if (strncmp((char*)camera_id, FD_DEVICE_PREFIX, strlen(FD_DEVICE_PREFIX)) == 0) {
			found_fd = true;
			break;
		}
	}

	if (!found_fd) {
		printf("Failed to find filter drive interface for the camera.\n");
		exit(1);
	}
	printf("Camera found: %s\n", camera_id);

	//attempt to open camera
	retval = J_Camera_Open(*factory, camera_id, camera);
	if (retval != J_ST_SUCCESS) {
		printf("Failed to open camera.\n");
		exit(1);
	}
}

image_collector::image_collector() : num_frames(0)
{
}

void image_collector::stream_callback(J_tIMAGE_INFO *raw_image)
{
	//convert raw image
	J_tIMAGE_INFO converted_image;
	J_Image_Malloc(raw_image, &converted_image);
	J_Image_FromRawToImage(raw_image, &converted_image);

	//add to list
	this->converted_frames.push_back(converted_image);

	num_frames++;
}

void uv_close_all(FACTORY_HANDLE *factory, CAM_HANDLE *camera)
{
	J_Camera_SetValueString(*camera, PROP_TRIGGER_MODE, TRIGGER_MODE_OFF);
	J_Camera_Close(*camera);
	J_Factory_Close(*factory);
}

NODE_HANDLE uv_camera_exposure_time_property_node(CAM_HANDLE *camera)
{
	NODE_HANDLE exposure_time;
	J_Camera_GetNodeByName(*camera, NODE_EXPOSURE_TIME, &exposure_time);
	return exposure_time;
}

void uv_camera_exposure_time_minmax(CAM_HANDLE *camera, double *min, double *max)
{
	NODE_HANDLE exposure_node = uv_camera_exposure_time_property_node(camera);
	J_Node_GetMinDouble(exposure_node, min);
	J_Node_GetMaxDouble(exposure_node, max);
}

void uv_camera_raw_exposure_time_minmax(CAM_HANDLE *camera, long *min, long *max)
{
	//get node
	NODE_HANDLE node;
	J_Camera_GetNodeByName(*camera, NODE_RAW_EXPOSURE_TIME, &node);

	//get ranges
	int64_t int64val;
	J_Node_GetMinInt64(node, &int64val);
	*min = int64val;
	J_Node_GetMaxInt64(node, &int64val);
	*max = int64val;
}

double uv_camera_exposure_time(CAM_HANDLE *camera)
{
	//use the node interface for parameter getting. Attempting to get exposure time directly by name does not work properly (probably not meant to, multiple parameters belong to a single category)
	NODE_HANDLE exposure_node = uv_camera_exposure_time_property_node(camera);
	double exposure_time;
	J_Node_GetValueDouble(exposure_node, true, &exposure_time);
	return exposure_time;
}

long uv_camera_raw_exposure_time(CAM_HANDLE *camera)
{
	NODE_HANDLE node;
	J_Camera_GetNodeByName(*camera, NODE_RAW_EXPOSURE_TIME, &node);

	int64_t int64val;
	J_Node_GetValueInt64(node, true, &int64val);
	long exposure_time_raw = int64val;
	return exposure_time_raw;
}

void uv_camera_set_exposure_time(CAM_HANDLE *camera, double exposure_time)
{
	NODE_HANDLE exposure_node = uv_camera_exposure_time_property_node(camera);
	J_Node_SetValueDouble(exposure_node, true, exposure_time);
}

NODE_HANDLE uv_camera_gain_property_node(CAM_HANDLE *camera)
{
	NODE_HANDLE gain;
	J_Camera_GetNodeByName(*camera, NODE_GAIN, &gain);
	return gain;
}

void uv_camera_gain_minmax(CAM_HANDLE *camera, long *min, long *max)
{
	NODE_HANDLE gain_node = uv_camera_gain_property_node(camera);

	int64_t int64_val;

	J_Node_GetMinInt64(gain_node, &int64_val);
	*min = int64_val;

	J_Node_GetMaxInt64(gain_node, &int64_val);
	*max = int64_val;
}

long uv_camera_gain(CAM_HANDLE *camera)
{
	NODE_HANDLE gain_node = uv_camera_gain_property_node(camera);
	int64_t gain;
	J_Node_GetValueInt64(gain_node, true, &gain);
	return (long)gain;
}

void uv_camera_set_gain(CAM_HANDLE *camera, long input_gain)
{
	NODE_HANDLE gain_node = uv_camera_gain_property_node(camera);
	int64_t gain = (int64_t)input_gain;
	J_Node_SetValueInt64(gain_node, true, gain);
}

void uv_camera_set_pixel_format(CAM_HANDLE *camera, int8_t *pixel_format_string)
{
	NODE_HANDLE pixel_format_node;
	J_Camera_GetNodeByName(*camera, NODE_PIXEL_FORMAT, &pixel_format_node);
	J_Node_SetValueString(pixel_format_node, true, pixel_format_string);
}

void uv_camera_set_mono8bit(CAM_HANDLE *camera)
{
	uv_camera_set_pixel_format(camera, PIXELFORMAT_MONO_8);
}

void uv_camera_set_mono10bit(CAM_HANDLE *camera)
{
	uv_camera_set_pixel_format(camera, PIXELFORMAT_MONO_10_PACKED);
}

std::string uv_camera_pixel_format(CAM_HANDLE *camera)
{
	NODE_HANDLE pixel_format_node;
	J_Camera_GetNodeByName(*camera, NODE_PIXEL_FORMAT, &pixel_format_node);

	uint32_t size = 512;
	char pixel_format[512] = {0};
	J_Node_GetValueString(pixel_format_node, true, (int8_t*)pixel_format, &size);
	return std::string(pixel_format);
}

void uv_camera_histogram(J_tIMAGE_INFO *image, histogram_mono8_t *histogram)
{
	//clear histogram
	for (int i=0; i < HISTOGRAM_LENGTH; i++) {
		histogram->counts[i] = 0;
	}

	//count pixels for histogram
	int width = image->iSizeX;
	int height = image->iSizeY;
	PIXELVALUE pixel;
	for (int i=0; i < width; i++) {
		for (int j=0; j < height; j++) {
			POINT point;
			point.y = i;
			point.x = j;
			J_Image_GetPixel(image, &point, &pixel);

			//assume mono8bit
			int value = pixel.PixelValueUnion.Mono8Type.Value;
			histogram->counts[value]++;
		}
	}
}

std::vector<J_tIMAGE_INFO> uv_acquire_images(CAM_HANDLE *camera, size_t num_frames)
{
	SIZE frame_size = uv_frame_size(camera);
	THRD_HANDLE thread;

	//set up an image_collector instance as the callback object, with its stream_callback function as the callback function
	image_collector collector;
	J_IMG_CALLBACK_OBJECT callback_object = reinterpret_cast<J_IMG_CALLBACK_OBJECT>(&collector);
	J_IMG_CALLBACK_FUNCTION callback_function = reinterpret_cast<J_IMG_CALLBACK_FUNCTION>(&image_collector::stream_callback);

	//open stream, use collector's stream_callback() function for receivings the frames in the stream
	J_Image_OpenStream(*camera, 0, callback_object, callback_function, &thread);
	uv_camera_start_acquisition(camera);

	//wait until expected number of frames have been collected
	while (collector.get_num_frames() < num_frames) {
		//camera has adjustable exposure time, but the frameperiod seems to be the same, so sleeping for 100 ms is probably a safe bet
		Sleep(100);
	}

	uv_camera_stop_acquisition(camera);
	J_Image_CloseStream(thread);

	//get frames
	std::vector<J_tIMAGE_INFO> images = collector.get_images();

	//remove superfluous frames from the beginning of the array (keep only the last num_frames frames)
	size_t extra_frames = images.size() - num_frames;
	if (extra_frames > 0) {
		images.erase(images.begin() + extra_frames);
	}
	return images;
}

void uv_save_image(J_tIMAGE_INFO *image, std::string filename)
{
	J_Image_SaveFile(image, filename.c_str());
}

void uv_save_images(std::vector<J_tIMAGE_INFO> images, std::string filename, verbosity_t verbosity)
{
	if (verbosity > 0) printf("Saving files:\n");

	for (size_t i=0; i < images.size(); i++) {
		std::ostringstream frame_filename;
		frame_filename << filename;
		if (images.size() > 1) {
			//add frame number with one leading zero
			frame_filename << "_frame_" << std::setw(2) << std::setfill('0') << i;
		}
		frame_filename << ".tiff";
		if (verbosity > 0) printf("* %s\n", frame_filename.str().c_str());
		uv_save_image(&images[i], frame_filename.str());
	}
}

std::string uv_camera_add_timestamps_and_info(CAM_HANDLE *camera, std::string base_string)
{
	//get string containing camera parameters
	long gain = uv_camera_gain(camera);
	double exposure_time = uv_camera_exposure_time(camera);
	std::string pixel_format = uv_camera_pixel_format(camera);
	std::ostringstream properties;
	properties << "gain_" << gain << "_exposure_" << long(exposure_time) << "us_" << pixel_format;

	//get date string
	time_t curr_time = time(NULL);
	struct tm timeinfo;
	localtime_s(&timeinfo, &curr_time);
	const int TIME_LENGTH = 80;
	char time_str[TIME_LENGTH];
	strftime(time_str, TIME_LENGTH, "%Y-%m-%dT%H%M%S", &timeinfo);

	//hard-coded device string
	const std::string DEVICE_STRING = "CM-140GE-UV";

	return base_string + "_" + DEVICE_STRING + "_" + properties.str() + "_" + std::string(time_str);
}

void show_help(const char *program_name)
{
	//display initial description
	printf("\nUsage:\n");
	printf("%s [options] [output_filename_prefix]\n\n", program_name);
	printf("When acquiring images to file, the filename will be postfixed by a timestamp and various camera parameters.\n\n");
	printf("If no filename is specified, no images will be aquired, unless live preview option is set (images will be acquired, but shown live on the screen instead of being saved to file).\n\n");
	printf("Camera parameters will be printed to the terminal during no-acquisition mode.\n\n");
	printf("Default pixel format is monochrome 8bit (alternative would be 10 bit packed).\n\n");

	//options
	printf("Options:\n");
	printf("   --live-view Show live view of the camera stream, along with sliders\n"
	       "\t\tfor adjusting gain and exposure time. Can also\n"
	       "\t\twrite images to file.\n");
	printf("   --no-acquisition No images will be acquired and saved to file. Will also\n"
	       "\t\toverride live-view.\n");
	printf("   --show-parameters Equivalent with the above.\n");
	printf("   --num-frames=NUM_FRAMES Specify number of frames to acquire in image\n"
	       "\t\tacquisition model. Images will be saved to file with\n"
	       "\t\t_01, _02, ... postfix on the filename. Default is a single\n"
	       "\t\tframe.\n");
	printf("   --gain=GAIN Specify raw gain value of the camera. Is set to %ld if not\n"
	       "\t\tspecified. Gain ranges are specified by the camera, run the\n"
	       "\t\tprogram without any arguments to see the gain range.\n", DEFAULT_GAIN);
	printf("   --exposure-time=EXPOSURE_TIME Specify exposure time in microseconds.\n"
	       "\t\tIs set to %f if not specified. Exposure time ranges are\n"
	       "\t\tspecified by the camera, run the program without any arguments\n"
	       "\t\tto see the gain range.\n", DEFAULT_EXPOSURE_TIME);
	printf("-h,--help Show help\n");
}

void uv_camera_start_acquisition(CAM_HANDLE *camera)
{
	J_Camera_ExecuteCommand(*camera, CMD_ACQUISITION_START);
}

void uv_camera_stop_acquisition(CAM_HANDLE *camera)
{
	J_Camera_ExecuteCommand(*camera, CMD_ACQUISITION_STOP);
}

void copy_image(J_tIMAGE_INFO *dst, J_tIMAGE_INFO *src)
{
	//set variables
	dst->iPixelType = src->iPixelType;
	dst->iSizeX = src->iSizeX;
	dst->iSizeY = src->iSizeY;
	dst->iImageSize = src->iImageSize;
	dst->iTimeStamp = src->iTimeStamp;
	dst->iMissingPackets = src->iMissingPackets;
	dst->iAnnouncedBuffers = src->iAnnouncedBuffers;
	dst->iQueuedBuffers = src->iQueuedBuffers;
	dst->iOffsetX = src->iOffsetX;
	dst->iOffsetY = src->iOffsetY;
	dst->iAwaitDelivery = src->iAwaitDelivery;
	dst->iBlockId = src->iBlockId;
	dst->iPaddingX = src->iPaddingX;
	dst->iImageStatus = src->iImageStatus;

	//copy data
	memcpy(dst->pImageBuffer, src->pImageBuffer, src->iImageSize);
}
