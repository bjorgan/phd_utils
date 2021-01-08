#ifndef CAMERA_CONTROL_WIDGET_H_DEFINED
#define CAMERA_CONTROL_WIDGET_H_DEFINED

#include <QWidget>
#include "uv_camera.h"

class QCloseEvent;
class QLabel;
class histogram_chart;
class QLineEdit;

/**
 * Callback object for displaying live stream from camera.
 **/
class liveview_callback {
	public:
		/**
		 * Constructor.
		 *
		 * \param liveview_window Window for displaying live view
		 * \param histogram_rect ROI for histogram calculation
		 **/
		liveview_callback(VIEW_HANDLE *liveview_window, RECT histogram_rect);

		/**
		 * Callback function for receiving raw images from stream and displaying them in the live view window.
		 *
		 * \param raw_frame Raw frame from stream
		 **/
		void callback_function(J_tIMAGE_INFO *raw_image);

		/**
		 * Get current frame.
		 *
		 * \return Copy of current frame. Will have to be deallocated elsewhere.
		 **/
		J_tIMAGE_INFO get_current_frame();
	private:
		///Live window
		VIEW_HANDLE *liveview_window;
		///Container for converted images
		J_tIMAGE_INFO converted_image;
		///Whether has been allocated
		bool converted_image_allocated;
		///Mutex-like thing for ensuring that extraction of current frame will happen safely
		CRITICAL_SECTION converted_image_mutex;
};

/**
 * Starts a live view against the UV camera, and displays widgets for controlling the camera. Shuts down the live view and closes the factory and camera when exited.
 **/
class camera_control : public QWidget {
	Q_OBJECT
	public:
		/**
		 * Constructor. Creates the widget layout and starts the live preview
		 * from the camera.
		 *
		 * Layout contains:
		 *  - Slider for controlling exposure time, and associated
		 *  labels for min/max and current value
		 *  - The same for controlling raw gain
		 *  - Histogram, which is updated once every second unless the sliders are changed (will then change on property change)
		 *
		 * \param factory JAI SDK factory, created elsewhere
		 * \param camera Camera handle, created elsewhere
		 * \param parent Parent widget
		 **/
		camera_control(FACTORY_HANDLE *factory, CAM_HANDLE *camera, QWidget *parent = NULL);
	public slots:
		/**
		 * Stop liveview.
		 **/
		void stop_liveview();
		/**
		 * Display a live view of the camera stream. Will call exit(1)
		 * on failure.  Starts the liveview in a separate thread, and
		 * exits the function. Liveview has to be stopped using
		 * uv_stop_live_view(). The window is displayed using weird
		 * windows magic, not Qt.
		 *
		 * liveview_callback class above is set up as a callback
		 * function on the stream that is being started.
		 **/
		void start_liveview();
		/**
		 * Set exposure time of camera. Raw exposure time will be
		 * converted to us and issued to camera.
		 *
		 * \param raw_exposure_time Raw exposure time (arbitrary units)
		 **/
		void set_raw_exposure_time(int raw_exposure_time);
		/**
		 * Set raw gain.
		 *
		 * \param gain Gain value
		 **/
		void set_gain(int gain);

		/**
		 * Update displayed histogram with current values from
		 * liveview. Obtains the current frame, updates histogram
		 * values.
		 **/
		void update_histogram_from_liveview();

		/**
		 * Get current frame and save it to disk. Uses the string in
		 * filename_input as prefix filename, and the working directory
		 * as directory.
		 **/
		void save_current_frame();
	signals:
		/**
		 * Emitted whenever gain or exposure times are changed, from set_raw_exposure_time() and set_gain().
		 **/
		void camera_properties_changed();
	protected:
		/**
		 * Reimplemented from QWidget::closeEvent(), stops the liveview and destroys factory and camera objects
		 * before closing the window.
		 **/
		void closeEvent(QCloseEvent *event);
	private:
		///JAI SDK factory handle
		FACTORY_HANDLE *factory;
		///Camera handle
		CAM_HANDLE *camera;
		///Absolute exposure time ranges (us)
		double min_exposure, max_exposure;
		///Raw exposure time ranges
		long min_raw_exposure, max_raw_exposure;
		///Raw gain ranges
		long min_gain, max_gain;
		///Displays current exposure time
		QLabel *exposure_time_label;
		///Displays current gain
		QLabel *gain_label;
		///Window displaying live view
		VIEW_HANDLE liveview_window;
		///Thread handling the live view stream
		THRD_HANDLE liveview_thread;
		///Whether live view is running
		bool liveview_running;
		///Callback object for stream
		liveview_callback *stream_callback;
		///Histogram values
		histogram_mono8_t histogram;
		///Histogram chart
		histogram_chart *histogram_displayer;
		///Filename input
		QLineEdit *filename_input;
};

#endif
