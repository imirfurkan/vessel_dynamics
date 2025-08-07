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
find_package(Eigen3 REQUIRED)

# 2) Linting (skip copyright checks for now)
if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  set(ament_cmake_copyright_FOUND TRUE)
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

# 3) Declare the C++ node executable
add_executable(vessel_kinematics_node src/vessel_kinematics_node.cpp)

# 4) Link that executable against the ROS client libraries
ament_target_dependencies(vessel_kinematics_node
  rclcpp
  geometry_msgs
  nav_msgs
  tf2_ros
  tf2
)

# 5) Install the node so it shows up in 
#    install/<pkg-name>/lib/<pkg-name>/
install(
  TARGETS vessel_kinematics_node
  RUNTIME DESTINATION lib/${PROJECT_NAME}
)

# 6) Install our data files (URDF, launch, rviz) into
#    install/<pkg-name>/share/<pkg-n  ame>/
install(
  DIRECTORY
    launch
    rviz
    urdf
  DESTINATION share/${PROJECT_NAME}
)


# for vessel dynamics node
find_package(Eigen3 REQUIRED)

add_executable(vessel_dynamics_node src/vessel_dynamics_node.cpp)

ament_target_dependencies(vessel_dynamics_node
  rclcpp
  geometry_msgs
  nav_msgs
  tf2
  tf2_ros
)
target_link_libraries(vessel_dynamics_node Eigen3::Eigen)

install(TARGETS vessel_dynamics_node
  DESTINATION lib/${PROJECT_NAME})

# milliampere dynamics node

add_executable(milliampere_dynamics_node src/milliampere_dynamics_node.cpp)

ament_target_dependencies(milliampere_dynamics_node
  rclcpp
  geometry_msgs
  nav_msgs
  tf2
  tf2_ros
  std_msgs
)

target_link_libraries(milliampere_dynamics_node Eigen3::Eigen)

install(TARGETS milliampere_dynamics_node
  DESTINATION lib/${PROJECT_NAME}
)


# 7) Let ament know we’re done
ament_package()
