cmake_minimum_required(VERSION 2.6)
project(hyperread_example C CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)

add_library(hyperread SHARED readimage.cpp readimage_stdin.cpp)
target_link_libraries(hyperread boost_regex explain)
install(TARGETS hyperread DESTINATION lib)
install(FILES readimage.h readimage_stdin.h DESTINATION include/hyperspectral)

configure_file(readimage.h hyperspectral/readimage.h COPYONLY)
configure_file(readimage_stdin.h hyperspectral/readimage_stdin.h COPYONLY)


add_executable(hyperread_example examples/read_and_print_image.cpp)
target_link_libraries(hyperread_example hyperread)

add_executable(header_reader header_reader_main.cpp ../common/getopt_long_helpers.cpp)
target_link_libraries(header_reader hyperread)
install(TARGETS header_reader DESTINATION bin)

if (WITH_TESTING)
	enable_testing()
	add_subdirectory(tests)
endif()
