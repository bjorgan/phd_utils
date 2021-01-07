#most readimage tests
add_executable(readimage-t readimage-t.cpp)
target_link_libraries(readimage-t CppUTest boost_regex hyperread)
add_test(NAME readimage COMMAND readimage-t)
configure_file(data/test_image.img data/test_image.img COPYONLY)
configure_file(data/test_image.hdr data/test_image.hdr COPYONLY)
configure_file(data/test_image_multiline_header.hdr data/test_image_multiline_header.hdr COPYONLY)

#test hyperspectral piping between processes
add_executable(read_image_from_stdout read_image_from_stdout.cpp)
target_link_libraries(read_image_from_stdout hyperread)
add_executable(write_image_to_stdout write_image_to_stdout.cpp)
target_link_libraries(write_image_to_stdout hyperread)
add_test(NAME hyperspectral_piping COMMAND sh -c "$<TARGET_FILE:write_image_to_stdout> | $<TARGET_FILE:read_image_from_stdout>")

#test hyperspectral stdin reading using the stdin API thing
add_executable(hyperspectral-stdin-no-accumulation hyperspectral-stdin-no-accumulation.cpp)
add_executable(hyperspectral-stdin-accumulation hyperspectral-stdin-accumulation.cpp)
target_link_libraries(hyperspectral-stdin-no-accumulation hyperread)
target_link_libraries(hyperspectral-stdin-accumulation hyperread)
add_test(NAME hyperspectral_piping_no_accumulation COMMAND sh -c "$<TARGET_FILE:write_image_to_stdout> | $<TARGET_FILE:hyperspectral-stdin-no-accumulation>")
add_test(NAME hyperspectral_piping_accumulation COMMAND sh -c "$<TARGET_FILE:write_image_to_stdout> | $<TARGET_FILE:hyperspectral-stdin-accumulation>")

#test whether stdin api is able to start at an arbitrary line (line 13 in this case)
add_test(NAME hyperspectral_piping_no_accumulation_arbitrary_line_start COMMAND sh -c "$<TARGET_FILE:write_image_to_stdout> | $<TARGET_FILE:hyperspectral-stdin-no-accumulation> 13")
