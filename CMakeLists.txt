 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(TudatTesting)

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(BUILD_STYLE "standalone")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
  set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
  # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIB_ROOT)
  set(LIB_ROOT "${CODEROOT}/lib")
endif()
if(NOT BIN_ROOT)
  set(BIN_ROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BIN_ROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIB_ROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIB_ROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BIN_ROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Define the install targets to create a distribution.
if(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)
    set(TUDAT_BUNDLE_DISTRIBUTION_PATH "${CODEROOT}")
endif(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)

if(NOT TEMPLATE_APPLICATION_DISTRIBUTION_PATH)
    set(TEMPLATE_APPLICATION_DISTRIBUTION_PATH 
        "${TUDAT_BUNDLE_DISTRIBUTION_PATH}/tudatApplications/templateApplication")
endif(NOT TEMPLATE_APPLICATION_DISTRIBUTION_PATH)

# Install Template Application files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${TEMPLATE_APPLICATION_DISTRIBUTION_PATH}/TemplateApplication"
        PATTERN ".DS_STORE" EXCLUDE
        PATTERN "CMakeLists.txt.user" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN ".bzr" EXCLUDE
)

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.45.0
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Find CSPICE library on local system.
find_package(Spice)

# Include CSpice directories.
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
endif( )

option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
if(NOT USE_NRLMSISE00)
  message(STATUS "NRLMSISE-00 disabled!")
  add_definitions(-DUSE_NRLMSISE00=0)
else()
  message(STATUS "NRLMSISE-00 enabled!")
  add_definitions(-DUSE_NRLMSISE00=1)
  # Find USE_NRLMSISE00 library on local system.
  find_package(NRLMSISE00)

  # Include NRLMSISE00 directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
  endif( )
endif( )

if(USE_CSPICE)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES tudat_spice_interface cspice)
endif()

if(USE_NRLMSISE00)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nrlmsise00)
endif()

list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES "")

if(USE_SOFA)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES sofa)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_sofa_interface )
 list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES tudat_earth_orientation )
endif()

if(USE_CSPICE)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES cspice)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_spice_interface )
endif()

if(USE_JSON)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nlohmann_json)
endif()

if(USE_NRLMSISE00)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nrlmsise00)
endif()

#if(USE_SPARTA)
#    add_definitions(-DUSE_SPARTA=1)
#    list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_sparta_interface)
#else()
#    add_definitions(-DUSE_SPARTA=0)
#endif()

list(APPEND TUDAT_APPLICATION_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_ground_stations tudat_propagators
   tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
   tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
   tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders
   tudat_filters tudat_basic_mathematics tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})

list(APPEND TUDAT_APPLICATION_ESTIMATION_LIBRARIES tudat_simulation_setup tudat_observation_models tudat_ground_stations tudat_acceleration_partials
   tudat_observation_partials tudat_estimatable_parameters tudat_orbit_determination tudat_propagators
   tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
   tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
   tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_root_finders
   tudat_filters tudat_basic_mathematics tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})


## Add aerobraking test
#add_executable(application_aerobrakingTest "${SRCROOT}/aerobrakingTest.cpp")
#setup_executable_target(application_aerobrakingTest "${SRCROOT}")
#target_link_libraries(application_aerobrakingTest ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

## Add Mars test
#add_executable(application_marsTest "${SRCROOT}/marsTest.cpp")
#setup_executable_target(application_marsTest "${SRCROOT}")
#target_link_libraries(application_marsTest ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

## Add conversion test
#add_executable(application_conversionTest "${SRCROOT}/conversionTest.cpp")
#setup_executable_target(application_conversionTest "${SRCROOT}")
#target_link_libraries(application_conversionTest ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add propagators comparison
add_executable(application_propagatorsComparison "${SRCROOT}/propagatorsComparison.cpp")
setup_executable_target(application_propagatorsComparison "${SRCROOT}")
target_link_libraries(application_propagatorsComparison ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(application_propagatorsComparisonExtra "${SRCROOT}/propagatorsComparisonExtra.cpp")
#setup_executable_target(application_propagatorsComparisonExtra "${SRCROOT}")
#target_link_libraries(application_propagatorsComparisonExtra ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

## Add tabulated atmosphere test
#add_executable(application_atmosphereTest "${SRCROOT}/atmosphereTest.cpp")
#setup_executable_target(application_atmosphereTest "${SRCROOT}")
#target_link_libraries(application_atmosphereTest tudat_aerodynamics tudat_input_output tudat_interpolators
#    tudat_basic_mathematics tudat_basic_astrodynamics ${Boost_LIBRARIES} )

# Add SPARTA test
#if(USE_SPARTA)
#    add_executable(application_spartaTest "${SRCROOT}/spartaTest.cpp")
#    setup_executable_target(application_spartaTest "${SRCROOT}")
#    target_link_libraries(application_spartaTest tudat_aerodynamics tudat_input_output tudat_basic_mathematics
#        tudat_basic_astrodynamics tudat_sparta_interface ${Boost_LIBRARIES} )
#endif()

## Add rotational motion test
#add_executable(application_rotationTest "${SRCROOT}/rotationTest.cpp")
#setup_executable_target(application_rotationTest "${SRCROOT}")
#target_link_libraries(application_rotationTest ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add quaternion propagation test
#add_executable(application_quaternionTest "${SRCROOT}/quaternionTest.cpp")
#setup_executable_target(application_quaternionTest "${SRCROOT}")
#target_link_libraries(application_quaternionTest ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )
