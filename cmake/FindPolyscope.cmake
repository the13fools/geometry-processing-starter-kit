# - Try to find the Polyscope library
# Once done this will define
#
#  POLYSCOPE_FOUND - system has POLYSCOPE 
#  POLYSCOPE_INCLUDE_DIR - **the** POLYSCOPE include directory
if(POLYSCOPE_FOUND)
    return()
endif()

find_path(POLYSCOPE_INCLUDE_DIR polyscope/polyscope.h
    HINTS
        ENV POLYSCOPE_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/polyscope
        ${CMAKE_SOURCE_DIR}/../polyscope
        ${CMAKE_SOURCE_DIR}/../tools/polyscope
        ${CMAKE_SOURCE_DIR}/../../polyscope
        /usr
        /usr/local
        ${CMAKE_SOURCE_DIR}/deps/polyscope
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(POLYSCOPE 
    "\npolyscope not found"
    POLYSCOPE_INCLUDE_DIR)
mark_as_advanced(POLYSCOPE_INCLUDE_DIR)
# "${POLYSCOPE_INCLUDE_DIR}/../"
add_subdirectory("${POLYSCOPE_INCLUDE_DIR}/../" "polyscope")
#add_subdirectory("/Users/hcsong/Workspaces/SGI22/drawing-with-nonmanifold-minimal-surface/TraceFrameField/deps/polyscope" "polyscope")