
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(TLDDC)
cmake_minimum_required (VERSION 3.1)
cmake_policy(SET CMP0057 NEW)
# -----------------------------------------------------------------------------
# CPP11
# -----------------------------------------------------------------------------
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED TRUE)
# -----------------------------------------------------------------------------
# DGtal
# -----------------------------------------------------------------------------
FIND_PACKAGE(DGtal 1.4 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")

# -----------------------------------------------------------------------------
# CLI11
# -----------------------------------------------------------------------------
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )

# -----------------------------------------------------------------------------
# bin
# -----------------------------------------------------------------------------

ADD_EXECUTABLE(logRec MainLogReconstruction3D IOHelper
    CylindricalPoint CylindricalCoordinateSystem)
TARGET_LINK_LIBRARIES(logRec ${DGTAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
