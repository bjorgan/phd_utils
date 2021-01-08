add_executable(calibration_io-t calibration_io-t.cpp ../src/calibration_io.cpp ../src/calibration.cpp ../../common/progress_bar.cpp)
target_link_libraries(calibration_io-t CppUTest hyperread spectral opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
add_test(NAME calibration_io COMMAND calibration_io-t)
