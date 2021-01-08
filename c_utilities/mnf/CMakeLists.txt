cmake_minimum_required(VERSION 2.6)
project(mnf C CXX)

#BLAS library
if (USE_MKL_LIBRARIES)
	message("Using MKL libraries.")
	set(BLAS_LIBRARIES mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
else()
	find_package(BLAS REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../common)
set(COMMON_SRC ../common/getopt_long_helpers.cpp ../common/string_helpers.cpp ../common/progress_bar.cpp)

#mnf library
add_library(mnf SHARED src/mnf.cpp src/mnf_statistics.cpp)
target_link_libraries(mnf ${BLAS_LIBRARIES} lapacke hyperread)

#mnf executable
add_executable(mnf-bin src/main.cpp ${COMMON_SRC})
set_target_properties(mnf-bin PROPERTIES OUTPUT_NAME mnf) #set mnf-bin binary name to mnf
target_link_libraries(mnf-bin mnf ${BLAS_LIBRARIES} lapacke hyperread)
configure_file(include/mnf.h hyperspectral/mnf.h COPYONLY)

#mnf-lbl executable
add_executable(mnf-lbl-bin src/main_linebyline.cpp ${COMMON_SRC})
set_target_properties(mnf-lbl-bin PROPERTIES OUTPUT_NAME mnf-linebyline)
target_link_libraries(mnf-lbl-bin mnf ${BLAS_LIBRARIES} lapacke hyperread)

#mnf utility
add_executable(mnf-info utilities/main.cpp ${COMMON_SRC})
target_link_libraries(mnf-info mnf hyperread)

#install libraries, binaries and header files
install(TARGETS mnf DESTINATION lib)
install(TARGETS mnf-bin DESTINATION bin)
install(TARGETS mnf-info DESTINATION bin)
install(TARGETS mnf-lbl-bin DESTINATION bin)
install(FILES include/mnf.h DESTINATION include/hyperspectral)

if (WITH_TESTING)
	enable_testing()
	add_subdirectory(tests/)
endif()
