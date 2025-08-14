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

# 3) milliampere dynamics node
add_executable(milliampere_dynamics_node src/milliampere_dynamics_node.cpp src/thrust_allocator.cpp)

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

install(TARGETS milliampere_dynamics_node
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
