
add_executable(main_bin arora-ge-ntru.cpp)

target_compile_options(main_bin PRIVATE -Wall -Werror -O2)

target_link_libraries(main_bin arora-ge-ntru)

set_target_properties(main_bin
  PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  OUTPUT_NAME arora-ge-ntru
)

#add_executable(keygen_bin keygen.cpp)

#target_compile_options(keygen_bin PRIVATE -Wall -Werror -O2)

#target_link_libraries(keygen_bin arora-ge-ntru)

#set_target_properties(keygen_bin
#  PROPERTIES 
#  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
#  OUTPUT_NAME keygen-ntru
#)

#add_executable(system_bin system.cpp)

#target_compile_options(system_bin PRIVATE -Wall -Werror -O2)

#target_link_libraries(system_bin arora-ge-ntru)

#set_target_properties(system_bin
#  PROPERTIES 
#  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
#  OUTPUT_NAME system-ntru
#)