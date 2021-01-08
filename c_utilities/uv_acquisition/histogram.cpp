#include "histogram.h"
#include <QLabel>
#include <QGridLayout>
#include <QPaintEvent>
#include <QTimer>

histogram_chart::histogram_chart(QWidget *parent) : QWidget(parent)
{
	image_width = HISTOGRAM_LENGTH;
	image_height = 100;
	image_data = new uchar[image_width*image_height*3]();

	histogram = new QLabel;
	min_value = new QLabel;
	max_value = new QLabel;
	max_value->setAlignment(Qt::AlignRight);
	min_value->setAlignment(Qt::AlignLeft);

	QGridLayout *layout = new QGridLayout(this);
	layout->addWidget(histogram, 0, 0, 1, 2);
	layout->addWidget(min_value, 1, 0);
	layout->addWidget(max_value, 1, 1);
}

void histogram_chart::update_histogram(histogram_mono8_t *histogram)
{
	//reset histogram image
	memset(image_data, 0, image_width*image_height*3);

	int min_value = 256;
	int max_value = -1;

	//set pixels in image to histogram values
	for (int i=0; i < image_width; i++) {
		int height = histogram->counts[i];

		if (height > 0) {
			int curr_value = i;
			if (curr_value < min_value) min_value = curr_value;
			if (curr_value > max_value) max_value = curr_value;
		}

		if (height > image_height) height = image_height;

		for (int j=0; j < height; j++) {
			int y = image_height - 1 - j;

			image_data[y*image_width*3 + i*3] = 255;
			image_data[y*image_width*3 + i*3+1] = 255;
			image_data[y*image_width*3 + i*3+2] = 255;
		}
	}

	this->min_value->setNum(min_value);
	this->max_value->setNum(max_value);

	update();
}

void histogram_chart::paintEvent(QPaintEvent *evt)
{
	Q_UNUSED(evt);

	QImage image(image_data, image_width, image_height, image_width*3, QImage::Format_RGB888); 
	histogram->setPixmap(QPixmap::fromImage(image));
	histogram->adjustSize();
}
