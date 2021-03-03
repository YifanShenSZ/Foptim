# FindFoptim
# -------
#
# Finds Foptim
#
# This will define the following variables:
#
#   Foptim_FOUND        -- True if the system has Foptim
#   Foptim_INCLUDE_DIRS -- The include directories for Foptim
#   Foptim_LIBRARIES    -- Libraries to link against
#
# and the following imported targets:
#
#   Foptim

# Find Foptim root
# Assume we are in ${FoptimROOT}/share/cmake/Foptim/FoptimConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(FoptimROOT "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# include directory
set(Foptim_INCLUDE_DIRS ${FoptimROOT}/include)

# library
add_library(Foptim STATIC IMPORTED)
set(Foptim_LIBRARIES Foptim)

# dependencies
list(APPEND Foptim_LIBRARIES ifcore matmul mkl_intel_lp64 mkl_intel_thread mkl_core)

# import
find_library(Foptim_LIBRARY Foptim PATHS "${FoptimROOT}/lib")
set_target_properties(Foptim PROPERTIES
    IMPORTED_LOCATION "${Foptim_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Foptim_INCLUDE_DIRS}"
    CXX_STANDARD 14
)