/*
// *  Copyright (c) 2018, Nagoya University
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of Autoware nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "ray_tracing.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <signal.h>

#include <ros/ros.h>
#include <std_msgs/String.h>
#include <sensor_msgs/PointCloud2.h>
#include <velodyne_pointcloud/point_types.h>
#include <velodyne_pointcloud/rawdata.h>
#include <tf/transform_listener.h>
#include <tf/message_filter.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <pcl_conversions/pcl_conversions.h>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

const int WINDOW_SIZE = 1000; //800;
const int MAX_DISTANCE = 500;
const int WIDTH=4096;
const int DEPTH=4096;
const int HEIGHT=20;
const double CELLSIZE = 0.3;

const double MAXHEIGHT_ABS_VAL = 2;
const double MINHEIGHT_ABS_VAL = 5;
const double HMAX=MAXHEIGHT_ABS_VAL;
const double HMIN=-MINHEIGHT_ABS_VAL;

double min_height = HMIN;
double max_height = HMAX;
double grid_width = WIDTH;
double grid_depth = DEPTH;
double grid_height = HEIGHT;
double cell_size = CELLSIZE;
double max_distance = MAX_DISTANCE;

pcl::PointCloud<velodyne_pointcloud::PointXYZIR> map;
tf::TransformListener *tf_listener;
ros::Publisher velodyne_pub;

pcl::PointCloud<pcl::PointXYZI> prev_points;
pcl::PointCloud<pcl::PointXYZI> pcl_out;

ros::Time prev_time;
double xo,yo,zo;

using namespace ground_field;
template class RayTracing<float>;

RayTracing<float> ray_tracing;

int stop_program = 0;

void signalHandler(int signum)
{
   if (signum == SIGINT) {
     stop_program = 1;
   }
}

void add_laser_points()
{
	bool added_rays = false;
	Point3D<float> origin(xo, yo, zo);

	std::cout << "To process " << ray_tracing.raysCount() << " rays " << "(" << pcl_out.points.size() << " points)" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    //Add all the rays
	ray_tracing.clearRays();
	for (auto p : pcl_out.points) {
		ray_tracing.addRay(Ray<float>(origin, Point3D<float>(p.x, p.y, p.z)));
	}
//    auto start = std::chrono::high_resolution_clock::now();
    //Perform parallel ray tracing
	#pragma omp parallel for default(shared)
    for (int i=0; i < ray_tracing.raysCount(); i++) {
        ray_tracing.computeRayIntersectionsZFilter2DGridEFLA(i);
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;

    //std::cout << "grimap: " << ray_tracing.getGridMap() << std::endl;
    std::cout << "Processed " << ray_tracing.raysCount() << " rays" << "(" << pcl_out.points.size() << " points)" << std::endl;
    std::cout << "Generated gridmap of " << ray_tracing.getGridMap().rows() << " rows, " << ray_tracing.getGridMap().cols() << " cols, " << ray_tracing.getGridMap().size() * sizeof(ray_tracing.getGridMap()(0,0)) << " bytes" << std::endl;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    std::cout << "Computed in " <<(double)microseconds/1000000.0 << " seconds" <<  std::endl;
}

void velodyneCallback(const pcl::PointCloud<pcl::PointXYZI>::ConstPtr& msg)
{
  static int count=0;

  std_msgs::Header header;
  pcl_conversions::fromPCL(msg->header,header);
  //printf("%f %f\n",header.stamp.toSec(),ros::Time::now().toSec());
  count++;

  tf::StampedTransform transform;
  if(prev_points.points.size()>0){
    try{
      tf_listener->waitForTransform("world","ndt_frame", prev_time/*header.stamp*/, ros::Duration(1));
      tf_listener->lookupTransform("world","ndt_frame",prev_time,transform);
    }catch(tf::TransformException ex){
      printf("old\n");
      return;
    }

    xo=transform.getOrigin().x();
    yo=transform.getOrigin().y();
    zo=transform.getOrigin().z();
    float dist = sqrt(xo*xo+yo*yo+zo*zo);
    printf("tf:x,y,z-> %f,%f,%f >>> %f %d\n", xo,yo,zo, dist, (dist > (WINDOW_SIZE/4)));

    pcl_out.clear();

    for(int i=0;i<prev_points.points.size(); i++){
      tf::Point pt(prev_points[i].x,prev_points[i].y,prev_points[i].z);
      tf::Point pt_world = transform * pt;
      pcl::PointXYZI wp;
      double distance=pt.x()*pt.x()+pt.y()*pt.y()+pt.z()*pt.z();
      if(distance<3*3)continue;
      wp.x=pt_world.x();
      wp.y=pt_world.y();
      wp.z=pt_world.z();
      wp.intensity=prev_points[i].intensity;

      pcl_out.push_back(wp);
    }

    pcl_out.header=prev_points.header;
    pcl_out.header.frame_id="world";
    velodyne_pub.publish(pcl_out);
  }
  prev_points=*msg;
  prev_time =header.stamp;

  add_laser_points();
}

//!
//! \brief getParam Get parameter from node handle
//! \param nh The nodehandle
//! \param param_name Key string
//! \param default_value Default value if not found
//! \return The parameter value
//!
template <typename T>
T getParam(const ros::NodeHandle& nh, const std::string& param_name, T default_value)
{
  T value;
  if (nh.hasParam(param_name))
  {
    nh.getParam(param_name, value);
  }
  else
  {
    ROS_WARN_STREAM("Parameter '" << param_name << "' not found, defaults to '" << default_value << "'");
    value = default_value;
  }
  return value;
}

void idle()
{
  while(!stop_program) {
	  usleep(100000);
	  ros::spinOnce();
  }
  throw NULL;
}

int main(int argc, char* argv[])
{
	  ros::init(argc, argv, "raytracing_velodyne");
	  ros::NodeHandle nh;

	  ros::Subscriber sub = nh.subscribe("velodyne_points", 10, velodyneCallback);
	  velodyne_pub = nh.advertise<pcl::PointCloud<pcl::PointXYZI > >("velodyne_points/world", 1);
	  tf_listener    = new tf::TransformListener();
	  sleep(2);

	  min_height = getParam(nh, "min_height", HMIN);
	  max_height = getParam(nh, "max_height", HMAX);
	  grid_width = getParam(nh, "grid_width", WIDTH);
	  grid_depth = getParam(nh, "grid_depth", DEPTH);
	  grid_height = getParam(nh, "grid_height", HEIGHT);
	  cell_size = getParam(nh, "cell_size", CELLSIZE);

	  printf("Running with min_height=%f, max_height=%f, grid_width=%f, grid_depth=%f, cell_size=%f\n", min_height, max_height, grid_width, grid_depth, cell_size);

	  //install the signal handler to be able to stop with CTRL-C
	  //signal(SIGINT, &signalHandler);

	  ray_tracing = RayTracing<float>(Point3D<float>(-grid_width/2,-grid_depth/2,-grid_height/2),
									  Point3D<float>(grid_width/2,grid_depth/2,grid_height/2),
									  Point3D<float>(cell_size,cell_size,cell_size));
	  //make the gridmap
	  ray_tracing.makeGridMapSize();

	  idle();
	  //ros::spin();

	  //ending
	  nh.deleteParam("min_height");
	  nh.deleteParam("max_height");
	  nh.deleteParam("grid_width");
	  nh.deleteParam("grid_depth");
	  nh.deleteParam("grid_height");
	  nh.deleteParam("cell_size");

	  return 0;
}
