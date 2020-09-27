CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
ADD_LIBRARY(
  density-tracer
  quaternion.cpp
)
INSTALL(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include/density-tracer
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
)
INSTALL(
  TARGETS density-tracer
  ARCHIVE
  DESTINATION lib
)
