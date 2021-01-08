cmake_minimum_required(VERSION 2.8.11)
project(uv-camera)
set(CMAKE_BUILD_TYPE Release)

#Qt setup
set(QT_MSVC_PATH "C:\\Qt\\Qt5.2.0\\5.2.0\\msvc2012_64_opengl")
set(CMAKE_PREFIX_PATH ${QT_MSVC_PATH})
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)

#JAI SDK path
set (JAI_SDK_PATH "C:/Program Files/JAI/SDK/library/CPP/")
include_directories(${JAI_SDK_PATH}/include)
link_directories(${JAI_SDK_PATH}/lib/Win64_x64)

#add getopt library
add_subdirectory(wingetopt)
include_directories(wingetopt/src)

#main exectuable
add_executable(uv_acquire main.cpp camera_control_widget.cpp histogram.cpp)
target_link_libraries(uv_acquire Jai_Factory wingetopt)
qt5_use_modules(uv_acquire Widgets Core)

#add_executable(test_gui main_test.cpp histogram.cpp)
#qt5_use_modules(test_gui Widgets Core)
