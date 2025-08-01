##
#  CMake script for time diffusion project
##

# Set the name of the project and target:
SET(TARGET "femffusion")


# Declare all source files the target consists of:
FILE(GLOB_RECURSE TARGET_SRC  "*.cc")
FILE(GLOB_RECURSE TARGET_INC  "*.h")
SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})

# Usually, you will not need to modify anything beyond this point...
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Are all dependencies fullfilled?
#
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_SLEPC)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_SLEPC = ON
which is required for this tutorial step."
    )
ENDIF()

project(femffusion CXX) 
# Enable OpenMP support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET (TARGET ${TARGET}.exe)
#include_directories(${CMAKE_SOURCE_DIR}/include)
DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
