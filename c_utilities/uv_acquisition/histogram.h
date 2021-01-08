#ifndef HISTOGRAM_H_DEFINED
#define HISTOGRAM_H_DEFINED

#include <QWidget>
#include "uv_camera.h"

class QLabel;
class QPaintEvent;

/**
 * Display histogram. Employs a QImage in representing the histogram values. See below for explanations.
 **/
class histogram_chart : public QWidget {
	Q_OBJECT
	public:
		/**
		 * Constructor. Creates histogram image and labels for displaying min and max.
		 *
		 * \param parent Parent widget
		 **/
		histogram_chart(QWidget *parent = NULL);
	public slots:
		/**
		 * Update image_data with new histogram data. For a histogram
		 * count of e.g. 67 for pixel value 10, image_data at column 10
		 * will be filled from pixel 0 to pixel 67. paintEvent is then
		 * implicitly called.
		 *
		 * \param histogram_data Histogram data
		 **/
		void update_histogram(histogram_mono8_t *histogram_data);
	protected:
		/**
		 * Reimplemented from QWidget. Constructs a QImage out of image_data,
		 * and displays it in the histogram QLabel.
		 *
		 * \param evt Event (unused)
		 **/
		void paintEvent(QPaintEvent *evt);
	private:
		///Width of histogram displayer image. Corresponds to the number of bins in the histogram.
		size_t image_width;
		///Height of histogram displayer image. Corresponds to the axis limit of the count, will cutoff at this value.
		size_t image_height;
		///Image data for histogram displayer image. Each column
		//corresponds to a pixel value in the histogram, and rows are
		//filled from the bottom and up to the histogram count for that
		//pixel value.
		uchar *image_data;
		///Displays the histogram displayer image.
		QLabel *histogram;
		///Displays current minimum and maximum pixel values
		QLabel *min_value, *max_value;
};



#endif
