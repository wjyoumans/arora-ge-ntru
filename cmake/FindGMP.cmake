find_path(GMP_INCLUDE_DIRS 
  NAMES gmp.h
  PATHS 
    $ENV{GMPDIR} 
    ${CMAKE_INSTALL_PREFIX}/include
  NO_DEFAULT_PATH
)

find_library(GMP_LIBRARIES
  NAMES gmp
  PATHS 
    $ENV{GMPDIR} 
    ${CMAKE_INSTALL_PREFIX}/lib 
  NO_DEFAULT_PATH
)

# If not found search default paths

find_path(GMP_INCLUDE_DIRS gmp.h)

find_library(GMP_LIBRARIES gmp)

# confirm GMP found

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP 
  DEFAULT_MSG
  GMP_INCLUDE_DIRS
  GMP_LIBRARIES
)

message(STATUS "Found GMP lib: ${GMP_LIBRARIES}")
