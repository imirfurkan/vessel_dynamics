cmake_minimum_required(VERSION 3.8)
project(vessel_kinematics)

# (Optional) extra warnings under GCC/Clang
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 1) Find ament macros & core ROS2 packages
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(tf2 REQUIRED)   # for tf2::Quaternion
find_package(std_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(osqp_vendor REQUIRED)


# 2) Linting (skip copyright checks for now)
if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  set(ament_cmake_copyright_FOUND TRUE)
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()


#################################################################
# --- Milliampere Dynamics Node ---
#################################################################

add_executable(milliampere_dynamics_node src/milliampere_dynamics_node.cpp)

ament_target_dependencies(milliampere_dynamics_node
  rclcpp
  geometry_msgs
  nav_msgs
  tf2
  tf2_ros
  std_msgs
  osqp_vendor
)

target_include_directories(milliampere_dynamics_node PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${osqp_vendor_INCLUDE_DIRS} # <-- ADD THIS LINE 

)

target_link_libraries(milliampere_dynamics_node Eigen3::Eigen osqp::osqp)

#################################################################
# --- Thrust Allocator Library & Node ---
#################################################################

add_library(vessel_lib src/thrust_allocator.cpp)
target_include_directories(vessel_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${osqp_vendor_INCLUDE_DIRS}
)
target_link_libraries(vessel_lib PUBLIC Eigen3::Eigen osqp::osqp)

add_executable(thrust_allocator_node src/thrust_allocator_node.cpp)

ament_target_dependencies(thrust_allocator_node
  rclcpp
  geometry_msgs
  std_msgs
)
target_link_libraries(thrust_allocator_node vessel_lib)


install(TARGETS
  milliampere_dynamics_node
  thrust_allocator_node
  DESTINATION lib/${PROJECT_NAME}
)

#################################################################
# --- Controller Node ---
#################################################################

add_executable(milliampere_controller_node src/milliampere_controller_node.cpp)

# CORRECTED: Only include your own package's headers here.
target_include_directories(milliampere_controller_node PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

ament_target_dependencies(milliampere_controller_node
  rclcpp
  geometry_msgs
  std_msgs
  nav_msgs
  tf2
  tf2_ros
)

# Link Eigen if you use it in this file
target_link_libraries(milliampere_controller_node
  Eigen3::Eigen
  vessel_lib            # if your node uses your library
)

install(TARGETS milliampere_controller_node
  DESTINATION lib/${PROJECT_NAME})

#################################################################
# --- Installation Rules ---
#################################################################

install(TARGETS
  milliampere_dynamics_node
  thrust_allocator_node
  milliampere_controller_node
  DESTINATION lib/${PROJECT_NAME}
)

# 4) Install non-source files (URDF, launch, rviz, etc.)
install(
  DIRECTORY
    launch
    rviz
    urdf
  DESTINATION share/${PROJECT_NAME}
)

# 5) Let ament know we’re done
ament_package()
