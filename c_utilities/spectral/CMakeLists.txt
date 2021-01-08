cmake_minimum_required(VERSION 2.6)
project(libspectral)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(spectral SHARED spectral.cpp)
install(TARGETS spectral DESTINATION lib)
install(FILES spectral.h DESTINATION include/spectral)

configure_file(spectral.h spectral/spectral.h COPYONLY)
