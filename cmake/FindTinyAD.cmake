# - Try to find the TinyAd library
# Once done this will define
#
#  TINYAD_FOUND - system has TINYAD 
#  TINYAD_INCLUDE_DIR - **the** TINYAD include directory
if(TINYAD_FOUND)
    return()
endif()

find_path(TINYAD_INCLUDE_DIR TinyAD/Scalar.hh
    HINTS
        ENV TINYAD_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/tinyad
        ${CMAKE_SOURCE_DIR}/../tinyad
        ${CMAKE_SOURCE_DIR}/../tools/tinyad
        ${CMAKE_SOURCE_DIR}/../../tinyad
        /usr
        /usr/local
        ${CMAKE_SOURCE_DIR}/deps/tinyad
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TinyAD 
    "\ntinyad not found"
    TINYAD_INCLUDE_DIR)
mark_as_advanced(TINYAD_INCLUDE_DIR)
# "${TINYAD_INCLUDE_DIR}/../"

if(TINYAD_INCLUDE_DIR)
    # Add subdirectory or do further processing
    add_subdirectory("${TINYAD_INCLUDE_DIR}/../" "tinyad")
else()
    message(FATAL_ERROR "TinyAd include directory not found. Please check your setup.")
endif()

# add_subdirectory("${TINYAD_INCLUDE_DIR}/../" "tinyad")
#add_subdirectory("/Users/hcsong/Workspaces/SGI22/drawing-with-nonmanifold-minimal-surface/TraceFrameField/deps/tinyad" "tinyad")