
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
FIND_PACKAGE(DGtal 1.3 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")

# -----------------------------------------------------------------------------
# Thread
# -----------------------------------------------------------------------------
FIND_PACKAGE(Threads)

# -----------------------------------------------------------------------------
# CLI11 et CSV
# -----------------------------------------------------------------------------
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )

# -----------------------------------------------------------------------------
# GSL
# -----------------------------------------------------------------------------

FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARIES})

# -----------------------------------------------------------------------------
# PCL
# -----------------------------------------------------------------------------

find_package(PCL 1.3 REQUIRED COMPONENTS common segmentation)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES})

# -----------------------------------------------------------------------------
# bin
# -----------------------------------------------------------------------------

ADD_EXECUTABLE(logRec MainLogReconstruction3D IOHelper
    CylindricalPoint CylindricalCoordinateSystem)
TARGET_LINK_LIBRARIES(logRec ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES})
