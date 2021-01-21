cmake_minimum_required(VERSION 2.8.3)
project(lidar_road_surface)

## Compile as C++14, supported in ROS Kinetic and newer
add_compile_options(-std=c++14 -march=native)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
    roscpp
    rospy
    std_msgs
    tf
    visualization_msgs
    message_generation
    pcl_ros
    sensor_msgs
    pcl_conversions
    velodyne_pointcloud
    grid_map_ros
    grid_map_cv
    grid_map_msgs
)

## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS filesystem date_time thread system)

find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Eigen3 QUIET)
if(NOT EIGEN3_FOUND)
  # Fallback to cmake_modules
  find_package(cmake_modules REQUIRED)
  find_package(Eigen REQUIRED)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})
  set(EIGEN3_LIBRARIES ${EIGEN_LIBRARIES})  # Not strictly necessary as Eigen is head only
  # Possibly map additional variables to the EIGEN3_ prefix.
  add_compile_options(-DUSE_EIGEN)
else()
  set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
  add_compile_options(-DUSE_EIGEN)
endif()

find_package(OpenGL REQUIRED QUIET)
find_package(GLUT REQUIRED)
find_package(GLU REQUIRED)
find_package(GLEW REQUIRED QUIET)

#find_package(PCL 1.7 REQUIRED COMPONENTS common features)
find_package(PCL 1.8 REQUIRED COMPONENTS gpu_octree common features)

find_package(Qt5Core)

find_package( OpenMP )
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

################################################
## Declare ROS messages, services and actions ##
################################################

###################################
## catkin specific configuration ##
###################################
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES   ${OPENGL_LIBRARIES} ${GLEW_LIBRARY}
  CATKIN_DEPENDS grid_map_msgs message_generation pcl_conversions pcl_ros roscpp rospy sensor_msgs std_msgs tf visualization_msgs
)


###########
## Build ##
###########
include_directories(include
  ${catkin_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIRS}
  ${GLEW_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  #nodes/ground_height_field/include
  nodes/ground_mapper
)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(CMAKE_CXX_FLAGS "-std=c++14 -O3 -g -Wall ${CMAKE_CXX_FLAGS}")

#add_executable(test_raytracing nodes/ground_height_field/test_raytracing.cpp)
#add_executable(test_raytracing_gridmap nodes/ground_height_field/test_raytracing_gridmap.cpp)
#add_executable(test_raytracing_velodyne nodes/ground_height_field/test_raytracing_velodyne.cpp)
#target_link_libraries(test_raytracing_velodyne ${catkin_LIBRARIES})

add_executable(ground_mapping
	nodes/ground_mapper/ground_mapper_node.cpp
	nodes/ground_mapper/mapping.cpp
	nodes/ground_mapper/grid.cpp)
target_link_libraries(ground_mapping ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${GLU_LIBRARY} ${GLEW_LIBRARY})
target_link_libraries(ground_mapping ${catkin_LIBRARIES} ${PCL_GPU_OCTREE_LIBRARIES} ${PCL_COMMON_LIBRARIES} ${PCL_FEATURES_LIBRARIES} ${Boost_LIBRARIES})

add_executable(pfm2pgm
	nodes/ground_mapper/pfm2pgm.cpp)
