cmake_minimum_required(VERSION 2.6)
project(hyper_picker C CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)
set(COMMON_SRC ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)

add_executable(hyper_picker main.cpp ${COMMON_SRC})
target_link_libraries(hyper_picker hyperread)
install(TARGETS hyper_picker DESTINATION bin)

add_executable(hyper_meaner main_hyper_meaner.cpp ${COMMON_SRC})
target_link_libraries(hyper_meaner hyperread)
install(TARGETS hyper_meaner DESTINATION bin)
