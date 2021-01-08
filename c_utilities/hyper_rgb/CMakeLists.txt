cmake_minimum_required(VERSION 2.6)
project(hyperread_example C CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(hyper_rgb main.cpp ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)
target_link_libraries(hyper_rgb hyperread opencv_highgui opencv_core opencv_imgproc opencv_imgcodecs)
install(TARGETS hyper_rgb DESTINATION bin)
