# 1) imports
from launch import LaunchDescription
from launch.substitutions import Command, FindExecutable, PathJoinSubstitution
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare

def generate_launch_description():
    # 2) locate our package on the filesystem
    pkg = FindPackageShare('vessel_kinematics')

    # 3) build the 'robot_description' parameter by running xacro on our URDF
    robot_desc = Command([
        # find the 'xacro' executable
        FindExecutable(name='xacro'), ' ',
        # path to vessel_prism.urdf.xacro inside install/share/...
        PathJoinSubstitution([pkg, 'urdf/milliampere.urdf.xacro'])
    ])

    # 4) return a LaunchDescription with the 3 nodes we want
    return LaunchDescription([

        # ————————————————————————————
        # A) robot_state_publisher
        #    reads the robot_description param and  
        #    broadcasts the static TFs for all links
        Node(
            package='robot_state_publisher',
            executable='robot_state_publisher',
            name='robot_state_publisher',
            parameters=[{'robot_description': robot_desc}],
            output='screen'
        ),
        
        # ————————————————————————————
        # B) RViz2
        #    launches rviz with the given .rviz config
        Node(
            package='rviz2',
            executable='rviz2',
            name='rviz2',
            output='screen',
            arguments=[
              '-d',
              PathJoinSubstitution([pkg, 'rviz/vessel.rviz'])
            ]
        ),

        # ————————————————————————————
        # C) milliampere_dynamics_node
        #    our simple kinematics, which publishes
        #    the dynamic map->base_link transform
        Node(
            package='vessel_kinematics',
            executable='milliampere_dynamics_node',
            name='milliampere_dynamics',
            output='screen'
        ),

        # ————————————————————————————
        # D) thrust_allocator_node
        #    
        #    
        Node(
            package='vessel_kinematics',
            executable='thrust_allocator_node',
            name='thrust_allocator',
            output='screen'
        ),
        
        # ————————————————————————————
        # E) pid_controller_node
        #    
        #    
        Node(
            package='vessel_kinematics',
            executable='pid_controller_node',
            name='pid_controller',
            output='screen'
        ),
        # # ————————————————————————————
        # # F) mpc_controller_node
        # #    
        # #    
        # Node(
        #     package='vessel_kinematics',
        #     executable='mpc_controller_node',
        #     name='mpc_controller',
        #     output='screen'
        # ),
        # # ————————————————————————————
        # G) disturbances_node
        #    
        #    
        Node(
            package='vessel_kinematics',
            executable='disturbances_node',
            name='disturbances',
            output='screen'
        ),
    ])
