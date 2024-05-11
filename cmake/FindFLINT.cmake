
# First search specified paths only

find_path(FLINT_INCLUDE_DIRS 
  NAMES flint.h
  PATHS 
    $ENV{FLINTDIR} 
    ${CMAKE_INSTALL_PREFIX}/include/flint
  NO_DEFAULT_PATH
)

find_library(FLINT_LIBRARIES
  NAMES flint
  PATHS 
    $ENV{FLINTDIR} 
    ${CMAKE_INSTALL_PREFIX}/lib
  NO_DEFAULT_PATH
)

# If not found search default paths

find_path(FLINT_INCLUDE_DIRS 
  NAMES flint.h
  PATH_SUFFIXES flint
)

find_library(FLINT_LIBRARIES NAMES flint)

# confirm FLINT found

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT
  DEFAULT_MSG
  FLINT_INCLUDE_DIRS
  FLINT_LIBRARIES
)

message(STATUS "Found FLINT lib: ${FLINT_LIBRARIES}")
