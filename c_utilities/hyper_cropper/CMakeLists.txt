cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)
set(COMMON_SRC ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)

add_executable(hyper_cropper main.cpp ${COMMON_SRC})
target_link_libraries(hyper_cropper hyperread)
install(TARGETS hyper_cropper DESTINATION bin)
