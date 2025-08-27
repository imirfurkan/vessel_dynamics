# generated from ament/cmake/core/templates/nameConfig.cmake.in

# prevent multiple inclusion
if(_vessel_kinematics_CONFIG_INCLUDED)
  # ensure to keep the found flag the same
  if(NOT DEFINED vessel_kinematics_FOUND)
    # explicitly set it to FALSE, otherwise CMake will set it to TRUE
    set(vessel_kinematics_FOUND FALSE)
  elseif(NOT vessel_kinematics_FOUND)
    # use separate condition to avoid uninitialized variable warning
    set(vessel_kinematics_FOUND FALSE)
  endif()
  return()
endif()
set(_vessel_kinematics_CONFIG_INCLUDED TRUE)

# output package information
if(NOT vessel_kinematics_FIND_QUIETLY)
  message(STATUS "Found vessel_kinematics: 0.0.0 (${vessel_kinematics_DIR})")
endif()

# warn when using a deprecated package
if(NOT "" STREQUAL "")
  set(_msg "Package 'vessel_kinematics' is deprecated")
  # append custom deprecation text if available
  if(NOT "" STREQUAL "TRUE")
    set(_msg "${_msg} ()")
  endif()
  # optionally quiet the deprecation message
  if(NOT vessel_kinematics_DEPRECATED_QUIET)
    message(DEPRECATION "${_msg}")
  endif()
endif()

# flag package as ament-based to distinguish it after being find_package()-ed
set(vessel_kinematics_FOUND_AMENT_PACKAGE TRUE)

# include all config extra files
set(_extras "")
foreach(_extra ${_extras})
  include("${vessel_kinematics_DIR}/${_extra}")
endforeach()
