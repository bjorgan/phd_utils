#include "camera_control_widget.h"
#include <QPushButton>
#include <QGridLayout>
#include <QCloseEvent>
#include <QSlider>
#include <QLabel>
#include <QTimer>
#include <QLineEdit>
#include "histogram.h"

camera_control::camera_control(FACTORY_HANDLE *factory, CAM_HANDLE *camera, QWidget *parent) : QWidget(parent), camera(camera), factory(factory)
{
	//create liveview
	this->liveview_running = false;
	start_liveview();

	//get raw and absolute exposure time ranges
	uv_camera_exposure_time_minmax(camera, &min_exposure, &max_exposure);
	uv_camera_raw_exposure_time_minmax(camera, &min_raw_exposure, &max_raw_exposure);

	//create exposure time slider
	exposure_time_label = new QLabel;
	QSlider *exposure_time = new QSlider(Qt::Horizontal);
	exposure_time->setMinimum(min_raw_exposure);
	exposure_time->setMaximum(max_raw_exposure);
	exposure_time->setValue(uv_camera_raw_exposure_time(camera));
	connect(exposure_time, SIGNAL(valueChanged(int)), SLOT(set_raw_exposure_time(int)));
	exposure_time_label->setText(QString::number(uv_camera_exposure_time(camera)));

	//get gain ranges
	uv_camera_gain_minmax(camera, &min_gain, &max_gain);

	//create gain slider
	gain_label = new QLabel;
	QSlider *gain = new QSlider(Qt::Horizontal);
	gain->setMinimum(min_gain);
	gain->setMaximum(max_gain);
	long current_gain = uv_camera_gain(camera);
	gain_label->setText(QString::number(uv_camera_gain(camera)));
	gain->setValue(current_gain);
	connect(gain, SIGNAL(valueChanged(int)), SLOT(set_gain(int)));

	//create filesave buttons
	filename_input = new QLineEdit;
	filename_input->setPlaceholderText("Filename prefix");
	QPushButton *save_button = new QPushButton("Save current frame to file");
	connect(save_button, SIGNAL(clicked()), SLOT(save_current_frame()));

	//histogram acquisition timer, update histogram once every second
	QTimer *timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), SLOT(update_histogram_from_liveview()));
	timer->start(1000);

	//also update histogram whenever gain or exposure times change
	connect(this, SIGNAL(camera_properties_changed()), SLOT(update_histogram_from_liveview()));

	histogram_displayer = new histogram_chart;

	//layout
	QGridLayout *layout = new QGridLayout(this);

	//add exposure time
	int row=0;
	layout->addWidget(new QLabel("Exposure time"), row, 0);
	layout->addWidget(exposure_time_label, row++, 1);
	layout->addWidget(exposure_time, row++, 0, 1, 3);

	//add exposure time min/max labels
	QLabel *max_exposure_label = new QLabel(QString::number(max_exposure));
	max_exposure_label->setAlignment(Qt::AlignRight);
	QLabel *min_exposure_label = new QLabel(QString::number(min_exposure));
	min_exposure_label->setAlignment(Qt::AlignLeft);
	layout->addWidget(min_exposure_label, row, 0);
	layout->addWidget(max_exposure_label, row++, 2);

	//add gain
	layout->addWidget(new QLabel("Raw gain"), row, 0);
	layout->addWidget(gain_label, row++, 1);
	layout->addWidget(gain, row++, 0, 1, 3);

	//add gain min/max labels
	QLabel *max_gain_label = new QLabel(QString::number(max_gain));
	max_gain_label->setAlignment(Qt::AlignRight);
	QLabel *min_gain_label = new QLabel(QString::number(min_gain));
	min_gain_label->setAlignment(Qt::AlignLeft);
	layout->addWidget(min_gain_label, row, 0);
	layout->addWidget(max_gain_label, row++, 2);

	//add histogram chart
	layout->addWidget(histogram_displayer, row++, 0, 1, 3);

	//add save widgets
	layout->addWidget(filename_input, row, 0, 1, 2);
	layout->addWidget(save_button, row, 2);
}

void camera_control::start_liveview()
{
	if (!this->liveview_running) {
		//get height and width of camera frames
		SIZE window_size = uv_frame_size(camera);

		//set window position
		POINT window_position = {0, 0};
		window_position.x = 10;
		window_position.y = 10;

		//open view window
		J_STATUS_TYPE retval = J_Image_OpenViewWindow("UV camera preview", &window_position, &window_size, &(liveview_window));
		if (retval != J_ST_SUCCESS) {
			printf("Failed to create preview window.\n");
			exit(1);
		}

		//ROI
		RECT rect = {0};
		rect.right = 100;
		rect.bottom = 100;

		//open stream and start acquisition
		stream_callback = new liveview_callback(&liveview_window, rect);
		J_IMG_CALLBACK_OBJECT callback_object = reinterpret_cast<J_IMG_CALLBACK_OBJECT>(stream_callback);
		J_IMG_CALLBACK_FUNCTION callback_function = reinterpret_cast<J_IMG_CALLBACK_FUNCTION>(&liveview_callback::callback_function);
		J_Image_OpenStream(*camera, 0, callback_object, callback_function, &liveview_thread);
		uv_camera_start_acquisition(camera);

		this->liveview_running = true;
	}
}

void camera_control::stop_liveview()
{
	if (this->liveview_running) {
		//stop acquisition
		uv_camera_stop_acquisition(camera);
		J_Image_CloseStream(liveview_thread);
		J_Image_CloseViewWindow(liveview_window);
		liveview_running = false;

		delete stream_callback;
	}
}

void camera_control::set_raw_exposure_time(int raw_exposure_time)
{
	long exposure_range = max_raw_exposure - min_raw_exposure;
	double fraction = (raw_exposure_time - min_raw_exposure)*1.0/((max_raw_exposure - min_raw_exposure)*1.0);

	double abs_exposure_time = fraction*(max_exposure - min_exposure) + min_exposure;
	exposure_time_label->setText(QString::number(abs_exposure_time));
	uv_camera_set_exposure_time(this->camera, abs_exposure_time);

	emit camera_properties_changed();
}

void camera_control::set_gain(int gain)
{
	uv_camera_set_gain(this->camera, gain);
	gain_label->setText(QString::number(gain));

	emit camera_properties_changed();
}

void camera_control::closeEvent(QCloseEvent *event)
{
	//stop liveview
	stop_liveview();

	//clean up after SDK
	uv_close_all(this->factory, this->camera);

	event->accept();
}

liveview_callback::liveview_callback(VIEW_HANDLE *liveview_window, RECT histogram_rect) : liveview_window(liveview_window), converted_image_allocated(false)
{
	InitializeCriticalSection(&converted_image_mutex);
}

void liveview_callback::callback_function(J_tIMAGE_INFO *raw_image)
{
	//show image
	J_Image_ShowImage(*liveview_window, raw_image);

	if (!converted_image_allocated) {
		J_Image_Malloc(raw_image, &converted_image);
		converted_image_allocated = true;
	}

	//convert raw to image
	EnterCriticalSection(&converted_image_mutex);
	J_Image_FromRawToImage(raw_image, &converted_image);
	LeaveCriticalSection(&converted_image_mutex);
}

J_tIMAGE_INFO liveview_callback::get_current_frame()
{
	EnterCriticalSection(&converted_image_mutex);

	//allocate ret frame
	J_tIMAGE_INFO ret_frame;
	J_Image_Malloc(&converted_image, &ret_frame);

	//copy image from converted frame to returned frame
	copy_image(&ret_frame, &converted_image);

	LeaveCriticalSection(&converted_image_mutex);
	return ret_frame;
}

void camera_control::update_histogram_from_liveview()
{
	if (liveview_running) {
		//get current frame
		J_tIMAGE_INFO current_frame = stream_callback->get_current_frame();

		//update histogram
		uv_camera_histogram(&current_frame, &histogram);

		//update histogram display
		histogram_displayer->update_histogram(&histogram);

		//free memory
		J_Image_Free(&current_frame);
	}
}

void camera_control::save_current_frame()
{
	if (liveview_running) {
		//get current frame
		std::vector<J_tIMAGE_INFO> frame;
		frame.push_back(stream_callback->get_current_frame());

		//add datetime postfix to filename
		std::string filename = filename_input->text().toStdString();
		filename = uv_camera_add_timestamps_and_info(camera, filename);

		//save to disk
		uv_save_images(frame, filename, VERBOSE_OUTPUT);

		//free memory
		J_Image_Free(&(frame[0]));
	}
}
