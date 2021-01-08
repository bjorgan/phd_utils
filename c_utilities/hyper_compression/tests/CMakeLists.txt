add_executable(mnf-compression-t mnf-compression-t.cpp)
target_link_libraries(mnf-compression-t CppUTest mnf-compression hyperread)
include_directories(..)
add_test(NAME mnf-compression-t COMMAND mnf-compression-t)
