add_executable(mnf-t mnf-t.cpp)
target_link_libraries(mnf-t CppUTest mnf)
include_directories(../src)
add_test(NAME mnf-t COMMAND mnf-t)
