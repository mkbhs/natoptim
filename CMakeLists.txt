cmake_minimum_required(VERSION 2.4.0)

project(NatOptim)

# Find source files
file(GLOB SOURCES src/*.cpp)

# Include header files
include_directories(include)

# Create shared library
add_library(${PROJECT_NAME} SHARED ${SOURCES})

# Install library
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

# Install library headers
file(GLOB HEADERS include/*.h)
install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})





#cmake_minimum_required(VERSION 3.5)
#project(natoptim VERSION 1.0.0)

#set(DEFAULT_BUILD_TYPE "Release")

#if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
#  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
#  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#endif()

#include(GNUInstallDirs)

#set(SOURCE_FILES src/randsearch.cpp)

#add_library(${PROJECT_NAME})

#target_include_directories(${PROJECT_NAME} PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)

#set_target_properties(${PROJECT_NAME} PROPERTIES
#    LINKER_LANGUAGE CXX 
#    VERSION ${PROJECT_VERSION}
#    SOVERSION 1)

#install(TARGETS ${PROJECT_NAME} EXPORT natoptimConfig
#    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})

#install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

#install(EXPORT natoptimConfig DESTINATION share/natoptim/cmake)

#export(TARGETS ${PROJECT_NAME} FILE natoptimConfig.cmake)

