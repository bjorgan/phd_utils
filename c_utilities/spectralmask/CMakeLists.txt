cmake_minimum_required(VERSION 2.6)
project(masking C CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)

set(COMMON_SRC ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)
add_library(masking SHARED src/masking.cpp)
add_executable(masking-bin src/main.cpp ${COMMON_SRC})
set_target_properties(masking-bin PROPERTIES OUTPUT_NAME hyper_mask)
target_link_libraries(masking-bin masking spectral hyperread)

install(TARGETS masking DESTINATION lib)
install(TARGETS masking-bin DESTINATION bin)
install(FILES src/masking.h DESTINATION include/hyperspectral)
