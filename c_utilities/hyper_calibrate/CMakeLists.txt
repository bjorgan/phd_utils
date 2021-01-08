cmake_minimum_required(VERSION 2.6)
project(hyper_calibrate)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/../common)

#remove the following two lines if this should be compiled without OpenCV support
add_definitions(-DWITH_OPENCV)
set(CMake_LIBS opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS regex)
include_directories(${Boost_INCLUDE_DIRS})

set(COMMON_SRC ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)
add_library(calibrate src/calibration.cpp src/calibration_io.cpp src/divide_by_max.cpp ${COMMON_SRC})
target_link_libraries(calibrate spectral hyperread)


#main calibration executable
add_executable(hyper_calibrate src/main.cpp)
target_link_libraries(hyper_calibrate ${CMake_LIBS} ${Boost_LIBRARIES} spectral hyperread calibrate)
install(TARGETS hyper_calibrate DESTINATION bin)

#spectral to calibration info conversion
add_executable(spectral_data_to_calibration_info utilities/spectral_data_to_calibration_info.cpp)
target_link_libraries(spectral_data_to_calibration_info ${CMake_LIBS} ${Boost_LIBRARIES} spectral hyperread calibrate)
set_target_properties(spectral_data_to_calibration_info PROPERTIES OUTPUT_NAME calibration_info_from_spectra)
install(TARGETS spectral_data_to_calibration_info DESTINATION bin)

#convenience scripts
configure_file(utilities/combine_spectra.sh combine_spectra COPYONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/combine_spectra DESTINATION bin/)

configure_file(utilities/duplicate_spectrum.sh duplicate_spectrum COPYONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/duplicate_spectrum DESTINATION bin/)

if (WITH_TESTING)
	enable_testing()
	add_subdirectory(tests)
endif()
