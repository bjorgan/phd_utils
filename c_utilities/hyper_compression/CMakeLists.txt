cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../common)

add_library(mnf-compression SHARED compression.cpp ../common/progress_bar.cpp)
target_link_libraries(mnf-compression mnf tar magic)

add_executable(mnf-compress main.cpp ../common/getopt_long_helpers.cpp)
target_link_libraries(mnf-compress mnf-compression hyperread)

install(TARGETS mnf-compression DESTINATION lib)
install(TARGETS mnf-compress DESTINATION bin)
install(FILES compression.h DESTINATION include/hyperspectral)

if (WITH_TESTING)
	enable_testing()
	add_subdirectory(tests/)
endif()
