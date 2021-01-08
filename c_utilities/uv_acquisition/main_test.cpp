#include <QApplication>
#include "histogram.h"

int main(int argc, char *argv[])
{
	QApplication qapp(argc, argv);
	histogram_chart histogram;
	histogram.show();
	return qapp.exec();
}
