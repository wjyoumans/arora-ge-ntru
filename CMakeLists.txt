cmake_minimum_required(VERSION 3.11)
project(arora-ge-ntru LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include_directories(./include)
#include_directories(/usr/local/include/flint)

add_subdirectory(src)
add_subdirectory(apps)
