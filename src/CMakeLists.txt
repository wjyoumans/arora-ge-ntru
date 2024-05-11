
find_package(GMP REQUIRED)
find_package(FLINT REQUIRED)

add_library(arora-ge-ntru SHARED "")

target_include_directories(arora-ge-ntru
  PUBLIC
  ${GMP_INCLUDE_DIRS}
  ${FLINT_INCLUDE_DIRS}
  #"${CMAKE_INSTALL_PREFIX}/include/flint"
  
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/include/arora-ge-ntru
)

target_sources(arora-ge-ntru
  PRIVATE
    keygen.cpp
    system.cpp
    recover.cpp
    extras.cpp
)

target_compile_options(arora-ge-ntru PRIVATE -Wall -Werror -O2)

target_link_libraries(arora-ge-ntru
  PUBLIC
    ${GMP_LIBRARIES}
    ${FLINT_LIBRARIES}
)

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_target_properties(arora-ge-ntru
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
)
#set_target_properties(arora-ge-ntru
#  PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
#)
